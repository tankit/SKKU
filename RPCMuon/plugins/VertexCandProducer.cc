#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <string>
#include <fstream>

class VertexCandProducer : public edm::EDFilter
{
public:
  VertexCandProducer(const edm::ParameterSet& pset);
  ~VertexCandProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const;
  inline double particleMass(const unsigned int pdgId) const;
  inline int signedPdgId(const unsigned int absPdgId, const int charge) const;

private:
  edm::InputTag trackLabel_;

  unsigned int pdgId_, leg1Id_, leg2Id_;
  double mass1_, mass2_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  edm::InputTag veto_muonLabel_;
  double veto_maxDeltaR_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

  unsigned int minNumber_, maxNumber_;

  //const TrackerGeometry* trackerGeom_;
  const MagneticField* bField_;
};

VertexCandProducer::VertexCandProducer(const edm::ParameterSet& pset)
{
  edm::ParameterSet vetoMuonPSet = pset.getParameter<edm::ParameterSet>("vetoMuon");
  veto_muonLabel_ = vetoMuonPSet.getParameter<edm::InputTag>("src");
  veto_maxDeltaR_ = vetoMuonPSet.getParameter<double>("maxDeltaR");

  edm::ParameterSet trackPSet = pset.getParameter<edm::ParameterSet>("track");
  trackLabel_ = trackPSet.getParameter<edm::InputTag>("src");
  cut_minPt_ = trackPSet.getParameter<double>("minPt");
  cut_maxEta_ = trackPSet.getParameter<double>("maxEta");
  cut_trackChi2_ = trackPSet.getParameter<double>("chi2");
  cut_trackNHit_  = trackPSet.getParameter<int>("nHit");
  cut_trackSignif_ = trackPSet.getParameter<double>("signif");
  cut_DCA_ = trackPSet.getParameter<double>("DCA");

  edm::ParameterSet vertexPSet = pset.getParameter<edm::ParameterSet>("vertex");
  cut_vertexChi2_ = vertexPSet.getParameter<double>("chi2");
  cut_minLxy_ = vertexPSet.getParameter<double>("minLxy");
  cut_maxLxy_ = vertexPSet.getParameter<double>("maxLxy");
  cut_vtxSignif_ = vertexPSet.getParameter<double>("signif");

  pdgId_ = pset.getParameter<unsigned int>("pdgId");
  leg1Id_ = pset.getParameter<unsigned int>("leg1Id");
  leg2Id_ = pset.getParameter<unsigned int>("leg2Id");
  rawMassMin_ = pset.getParameter<double>("rawMassMin");
  rawMassMax_ = pset.getParameter<double>("rawMassMax");
  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  mass1_ = particleMass(leg1Id_);
  mass2_ = particleMass(leg2Id_);

  produces<reco::VertexCompositeCandidateCollection>();
}

bool VertexCandProducer::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef reco::VertexCompositeCandidateCollection VCCandColl;

  std::auto_ptr<VCCandColl> decayCands(new VCCandColl);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel("offlineBeamSpot", beamSpotHandle);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByLabel(trackLabel_, trackHandle);

  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  bField_ = bFieldHandle.product();

  //edm::ESHandle<TrackerGeometry> trackerGeomHandle;
  //eventSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
  //trackerGeom_ = trackerGeomHandle.product();

  edm::ESHandle<GlobalTrackingGeometry> glbTkGeomHandle;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(glbTkGeomHandle);
  //glbTkGeom_ = glbTkGeomHandle.product();

  double vetoMuonEta = -999, vetoMuonPhi = 0;
  edm::Handle<reco::MuonCollection> vetoMuonHandle;
  if ( event.getByLabel(veto_muonLabel_, vetoMuonHandle) )
  {
    const reco::Muon& vetoMuon = vetoMuonHandle->at(0);
    vetoMuonEta = vetoMuon.eta();
    vetoMuonPhi = vetoMuon.phi();
  }

  for ( int i=0, n=trackHandle->size(); i<n; ++i )
  {
    TrackRef trackRef1(trackHandle, i);
    // Positive particle in the 1st index (pi+, proton, K+...)
    if ( trackRef1->charge() < 0 ) continue;
    if ( !isGoodTrack(trackRef1, beamSpotHandle.product()) ) continue;
    if ( vetoMuonEta != -999 and deltaR(vetoMuonEta, vetoMuonPhi, trackRef1->eta(), trackRef1->phi()) < veto_maxDeltaR_ ) continue;

    TransientTrack transTrack1(*trackRef1, bField_, glbTkGeomHandle);
    if ( !transTrack1.impactPointTSCP().isValid() ) continue;
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();

    for ( int j=0; j<n; ++j )
    {
      TrackRef trackRef2(trackHandle, j);
      // Negative particle in the 2nd index (pi-, anti-proton, K-...)
      if ( trackRef2->charge() > 0 ) continue;
      if ( !isGoodTrack(trackRef2, beamSpotHandle.product()) ) continue;
      if ( vetoMuonEta != -999 and deltaR(vetoMuonEta, vetoMuonPhi, trackRef2->eta(), trackRef2->phi()) < veto_maxDeltaR_ ) continue;

      TransientTrack transTrack2(*trackRef2, bField_, glbTkGeomHandle);
      if ( !transTrack2.impactPointTSCP().isValid() ) continue;
      FreeTrajectoryState ipState2 = transTrack2.impactPointTSCP().theState();

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(ipState1, ipState2);
      if ( !cApp.status() ) continue;
      const float dca = fabs(cApp.distance());
      if ( dca < 0. || dca > cut_DCA_ ) continue;
      GlobalPoint cxPt = cApp.crossingPoint();
      if (std::hypot(cxPt.x(), cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue; 

      TrajectoryStateClosestToPoint caState1 = transTrack1.trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint caState2 = transTrack2.trajectoryStateClosestToPoint(cxPt);
      if ( !caState1.isValid() or !caState2.isValid() ) continue;

      double mass1 = mass1_, mass2 = mass2_;
      int leg1Id = leg1Id_, leg2Id = leg2Id_;
      if ( leg1Id != leg2Id )
      {
        if ( caState1.momentum().mag() > caState2.momentum().mag() )
        {
          if ( mass1 < mass2 )
          {
            std::swap(mass1, mass2);
            std::swap(leg1Id, leg2Id);
          }
        }
        else
        {
          if ( mass1 > mass2 )
          {
            std::swap(mass1, mass2);
            std::swap(leg1Id, leg2Id);
          }
        }
      }
 
      const double rawEnergy = std::hypot(caState1.momentum().mag(), mass1) 
                             + std::hypot(caState2.momentum().mag(), mass2);
      const double rawMass = sqrt(rawEnergy*rawEnergy - (caState1.momentum()+caState2.momentum()).mag2());
      if ( rawMassMin_ > rawMass or rawMassMax_ < rawMass ) continue;

      // Build Vertex
      std::vector<TransientTrack> transTracks;
      transTracks.push_back(transTrack1);
      transTracks.push_back(transTrack2);
      KalmanVertexFitter fitter(true);
      TransientVertex transVertex = fitter.vertex(transTracks);

      if ( !transVertex.isValid() or transVertex.totalChiSquared() < 0. ) continue;

      const reco::Vertex vertex = transVertex;
      if ( vertex.normalizedChi2() > cut_vertexChi2_ ) continue;

      std::vector<TransientTrack> refittedTracks;
      if ( transVertex.hasRefittedTracks() ) refittedTracks = transVertex.refittedTracks();

      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;

      GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z());
      GlobalPoint beamSpotPos(beamSpotHandle->position().x(),
                              beamSpotHandle->position().y(),
                              beamSpotHandle->position().z());

      SMatrixSym3D totalCov = beamSpotHandle->rotatedCovariance3D() + vertex.covariance();
      SVector3 distanceVector(vertex.x() - beamSpotPos.x(), vertex.y() - beamSpotPos.y(), 0.);

      double rVtxMag = ROOT::Math::Mag(distanceVector);
      double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / rVtxMag;
      if( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ ) continue;

      // Cuts finished, now we create the candidates and push them back into the collections.
      
      std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
      std::auto_ptr<TrajectoryStateClosestToPoint> traj2;

      if ( refittedTracks.empty() )
      {
        traj1.reset(new TrajectoryStateClosestToPoint(transTrack1.trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(transTrack2.trajectoryStateClosestToPoint(vtxPos)));
      }
      else
      {
        TransientTrack* refTrack1 = 0, * refTrack2 = 0;
        for ( std::vector<TransientTrack>::iterator refTrack = refittedTracks.begin();
              refTrack != refittedTracks.end(); ++refTrack )
        {
          if ( refTrack->track().charge() > 0 ) refTrack1 = &*refTrack;
          else refTrack2 = &*refTrack;
        }
        if ( refTrack1 == 0 or refTrack2 == 0 ) continue;
        traj1.reset(new TrajectoryStateClosestToPoint(refTrack1->trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(refTrack2->trajectoryStateClosestToPoint(vtxPos)));
      }
      if( !traj1->isValid() or !traj2->isValid() ) continue;

      GlobalVector mom1(traj1->momentum());
      GlobalVector mom2(traj2->momentum());
      GlobalVector mom(mom1+mom2);

      //cleanup stuff we don't need anymore
      traj1.reset();
      traj2.reset();

      Particle::Point vtx(vertex.x(), vertex.y(), vertex.z());
      const Vertex::CovarianceMatrix vtxCov(vertex.covariance());
      double vtxChi2(vertex.chi2());
      double vtxNdof(vertex.ndof());

      const double candE1 = hypot(mom1.mag(), mass1);
      const double candE2 = hypot(mom2.mag(), mass2);

      Particle::LorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
      if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

      RecoChargedCandidate cand1(trackRef1->charge(), Particle::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1), vtx);
      RecoChargedCandidate cand2(trackRef2->charge(), Particle::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2), vtx);
      cand1.setTrack(trackRef1);
      cand2.setTrack(trackRef2);
      const int pdgId1 = signedPdgId(leg1Id, trackRef1->charge());
      const int pdgId2 = signedPdgId(leg2Id, trackRef2->charge());
      cand1.setPdgId(pdgId1);
      cand2.setPdgId(pdgId2);
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      cand->addDaughter(cand1);
      cand->addDaughter(cand2);

      cand->setPdgId(pdgId_);
      AddFourMomenta addP4;
      addP4.set(*cand);

      decayCands->push_back(*cand);
      
    }
  }

  const unsigned int nCands = decayCands->size();
  event.put(decayCands);

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}

bool VertexCandProducer::isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const
{
  const static reco::TrackBase::TrackQuality trackQual = reco::TrackBase::qualityByName("loose");
  if ( !track->quality(trackQual) ) return false;
  if ( track->normalizedChi2() >= cut_trackChi2_ ) return false;
  if ( track->numberOfValidHits() < cut_trackNHit_ ) return false;
  if ( track->pt() < cut_minPt_ or abs(track->eta()) > cut_maxEta_ ) return false;

  FreeTrajectoryState initialFTS = trajectoryStateTransform::initialFreeState(*track, bField_);
  TSCBLBuilderNoMaterial blsBuilder;
  TrajectoryStateClosestToBeamLine tscb( blsBuilder(initialFTS, *beamSpot) );
  if ( !tscb.isValid() ) return false;
  if ( tscb.transverseImpactParameter().significance() <= cut_trackSignif_ ) return false;

  return true;
}

double VertexCandProducer::particleMass(const unsigned int pdgId) const
{
  switch(pdgId)
  {
    case   11: return 0.0005109989; break;
    case   13: return 0.1056583715; break;
    case  211: return 0.13957018  ; break;
    case  321: return 0.493677    ; break;
    case 2212: return 0.938272013 ; break;
  }
  return 1e12; // Make this particle supermassive to fail mass range cut
} 

int VertexCandProducer::signedPdgId(const unsigned int pdgId, const int charge) const
{
  if ( pdgId == 11 or pdgId == 13 or pdgId == 15 ) return -pdgId*charge;
  return pdgId*charge;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexCandProducer);
