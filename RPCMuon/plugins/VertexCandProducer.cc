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

#include <string>
#include <fstream>

class VertexCandProducer : public edm::EDFilter
{
public:
  VertexCandProducer(const edm::ParameterSet& pset);
  ~VertexCandProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot);
  std::map<int, double> massMap_;
  
private:
  edm::InputTag trackLabel_;

  unsigned int pdgId_, leg1Id_, leg2Id_;
  double mass1_, mass2_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackIPSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minVtxDxy_, cut_minVtxSignif_;

  unsigned int minNumber_, maxNumber_;

  //const TrackerGeometry* trackerGeom_;
  const MagneticField* bField_;
};

VertexCandProducer::VertexCandProducer(const edm::ParameterSet& pset)
{
  edm::ParameterSet trackPSet = pset.getParameter<edm::ParameterSet>("track");
  trackLabel_ = trackPSet.getParameter<edm::InputTag>("src");
  cut_minPt_ = trackPSet.getParameter<double>("minPt");
  cut_maxEta_ = trackPSet.getParameter<double>("maxEta");
  cut_trackChi2_ = trackPSet.getParameter<double>("chi2");
  cut_trackNHit_  = trackPSet.getParameter<int>("nHit");
  cut_trackIPSignif_ = trackPSet.getParameter<double>("ipSignif");
  cut_DCA_ = trackPSet.getParameter<double>("DCA");

  edm::ParameterSet vertexPSet = pset.getParameter<edm::ParameterSet>("vertex");
  cut_vertexChi2_ = vertexPSet.getParameter<double>("chi2");
  cut_minVtxDxy_ = vertexPSet.getParameter<double>("dxy");
  cut_minVtxSignif_ = vertexPSet.getParameter<double>("vtxSignif");

  pdgId_ = pset.getParameter<unsigned int>("pdgId");
  leg1Id_ = pset.getParameter<unsigned int>("leg1Id");
  leg2Id_ = pset.getParameter<unsigned int>("leg2Id");
  rawMassMin_ = pset.getParameter<double>("rawMassMin");
  rawMassMax_ = pset.getParameter<double>("rawMassMax");
  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  massMap_[211] = 0.13957018;
  massMap_[2211] = 0.938272013;
  massMap_[321] = 0.497614;
  massMap_[3222] = 1.115683;

  mass1_ = massMap_[leg1Id_];
  mass2_ = massMap_[leg2Id_];

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

  for ( int i=0, n=trackHandle->size(); i<n; ++i )
  {
    TrackRef trackRef1(trackHandle, i);
    // Positive particle in the 1st index (pi+, proton, K+...)
    if ( trackRef1->charge() < 0 ) continue;
    if ( !isGoodTrack(trackRef1, beamSpotHandle.product()) ) continue;

    TransientTrack transTrack1(*trackRef1, bField_, glbTkGeomHandle);
    if ( !transTrack1.impactPointTSCP().isValid() ) continue;
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();

    for ( int j=0; j<n; ++j )
    {
      TrackRef trackRef2(trackHandle, j);
      // Negative particle in the 2nd index (pi-, anti-proton, K-...)
      if ( trackRef2->charge() > 0 ) continue;
      if ( !isGoodTrack(trackRef2, beamSpotHandle.product()) ) continue;

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

      const double rawEnergy = std::hypot(caState1.momentum().mag(), mass1_) 
                             + std::hypot(caState2.momentum().mag(), mass2_);
      const double rawMass = sqrt(rawEnergy*rawEnergy - (caState1.momentum()+caState2.momentum()).mag2());
      if ( rawMassMin_ > rawMass or rawMassMax_ < rawMass ) continue;

      // Build Vertex
      std::vector<TransientTrack> transTracks;
      transTracks.push_back(transTrack1);
      transTracks.push_back(transTrack2);
      KalmanVertexFitter fitter(false);
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
      if( rVtxMag < cut_minVtxDxy_ or rVtxMag / sigmaRvtxMag < cut_minVtxSignif_ ) continue;

      // Cuts finished, now we create the candidates and push them back into the collections.
      
      std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
      std::auto_ptr<TrajectoryStateClosestToPoint> traj2;

      traj1.reset(new TrajectoryStateClosestToPoint(transTrack1.trajectoryStateClosestToPoint(vtxPos)));
      traj2.reset(new TrajectoryStateClosestToPoint(transTrack2.trajectoryStateClosestToPoint(vtxPos)));
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

      const double candE1 = sqrt(mom1.mag2() + mass1_*mass1_);
      const double candE2 = sqrt(mom2.mag2() + mass2_*mass2_);

      Particle::LorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
      RecoChargedCandidate cand1(trackRef1->charge(), Particle::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1), vtx);
      RecoChargedCandidate cand2(trackRef2->charge(), Particle::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2), vtx);
      cand1.setTrack(trackRef1);
      cand2.setTrack(trackRef2);
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      cand->addDaughter(cand1);
      cand->addDaughter(cand2);

      if ( leg1Id_ != leg2Id_ and trackRef1->charge() < 0 ) cand->setPdgId(-pdgId_);
      else cand->setPdgId(pdgId_);

      decayCands->push_back(*cand);
      
    }
  }

  const unsigned int nCands = decayCands->size();
  event.put(decayCands);

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}

bool VertexCandProducer::isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot)
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
  if ( tscb.transverseImpactParameter().significance() <= cut_trackIPSignif_ ) return false;

  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexCandProducer);
