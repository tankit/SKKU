#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class FakeMuonAnalyzer : public edm::EDAnalyzer
{
public:
  FakeMuonAnalyzer(const edm::ParameterSet& pset);
  ~FakeMuonAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:

  const reco::Muon* findMatchedMuonByTrackRef(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle);
  const reco::Muon* findMatchedMuonByDR(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle);
  const reco::Muon* findMatchedMuonByDRDPt(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle);

  edm::InputTag muonLabel_;
  edm::InputTag vertexCandLabel_;
  std::vector<StringCutObjectSelector<reco::Muon, true>* > muonTypes_;
  StringCutObjectSelector<reco::VertexCompositeCandidate, true>* vertexCut_;
  double maxDR_, maxDPt_;

  TH1F* hEvent_, * hNCandOverlap_;
  TTree* tree_;
  int run_, lumi_, event_;
  int nMuon_, nVertexCand_;
  int legId_;
  int muonCharge_, trackCharge_;
  double vertexMass_, vertexPt_, vertexL3D_, vertexL2D_, vertexLxy_;
  double deltaR_, deltaPt_;
  math::XYZTLorentzVector muon_, track_;
  std::vector<int> muonTypeResults_;
};

FakeMuonAnalyzer::FakeMuonAnalyzer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getParameter<edm::InputTag>("muon");
  vertexCandLabel_ = pset.getParameter<edm::InputTag>("vertexCand");

  std::string vertexCut = pset.getParameter<std::string>("vertexCut");
  vertexCut_ = new StringCutObjectSelector<reco::VertexCompositeCandidate, true>(vertexCut);

  std::string matchBy = pset.getParameter<std::string>("match");
  if ( matchBy == "dR")
  {
    maxDR_ = pset.getParameter<double>("maxDR");
  }
  else if ( matchBy == "dRdPt" )
  {
    maxDR_ = pset.getParameter<double>("maxDR");
    maxDPt_ = pset.getParameter<double>("maxDPt");
  }

  edm::Service<TFileService> fs;
  hEvent_ = fs->make<TH1F>("hEvent", "Event count", 5, 0, 5);
  hEvent_->GetXaxis()->SetBinLabel(1, "Total");
  hEvent_->GetXaxis()->SetBinLabel(2, "Muon");
  hEvent_->GetXaxis()->SetBinLabel(3, "Kshort");
  hEvent_->GetXaxis()->SetBinLabel(4, "Overlap");

  hNCandOverlap_ = fs->make<TH1F>("hNCandOverlap", "Number of vertex candidate with track sharing", 10, 0, 10);

  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");

  tree_->Branch("nMuon", &nMuon_, "nMuon/I");
  tree_->Branch("nKshort", &nVertexCand_, "nKshort/I");

  tree_->Branch("vertexMass", &vertexMass_, "vertexMass/D");
  tree_->Branch("vertexPt", &vertexPt_, "vertexPt/D");
  tree_->Branch("vertexL3D", &vertexL3D_, "vertexL3D/D");
  tree_->Branch("vertexL2D", &vertexL2D_, "vertexL2D/D");
  tree_->Branch("vertexLxy", &vertexLxy_, "vertexLxy/D");

  tree_->Branch("deltaR" , &deltaR_ , "deltaR/D" );
  tree_->Branch("deltaPt", &deltaPt_, "deltaPt/D");
  tree_->Branch("legId", &legId_, "legId/I");
  tree_->Branch("muonCharge", &muonCharge_, "muonCharge/I");
  tree_->Branch("trackCharge", &trackCharge_, "trackCharge/I");
  tree_->Branch("muon" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &muon_ );
  tree_->Branch("track", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &track_);

  edm::ParameterSet muonTypes = pset.getParameter<edm::ParameterSet>("muonTypes");
  std::vector<std::string> muonTypeNames = muonTypes.getParameterNames();
  muonTypeResults_.resize(muonTypeNames.size());
  for ( int i=0, n=muonTypeNames.size(); i<n; ++i )
  {
    const std::string& cutName = muonTypeNames[i];
    std::string muonType = muonTypes.getParameter<std::string>(cutName);
    muonTypes_.push_back(new StringCutObjectSelector<reco::Muon, true>(muonType));
    tree_->Branch(Form("muonType_%s", cutName.c_str()), &(muonTypeResults_[0])+i, Form("muonType_%s/I", cutName.c_str()));
  }
}

void FakeMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  run_ = event.id().run();
  event_ = event.id().event();
  lumi_ = event.id().luminosityBlock();

  for ( int i=0, n=muonTypeResults_.size(); i<n; ++i ) muonTypeResults_[i] = -1;

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  edm::Handle<edm::View<reco::VertexCompositeCandidate> > vertexCandHandle;
  event.getByLabel(vertexCandLabel_, vertexCandHandle);

  nMuon_ = muonHandle->size();
  nVertexCand_ = vertexCandHandle->size();

  hEvent_->Fill(0);
  if ( nMuon_ != 0 ) hEvent_->Fill(1);
  if ( nVertexCand_ != 0 ) hEvent_->Fill(2);

  int nFakeMuon = 0;

  std::vector<int> overlappedCandIds;
  for ( int iVertexCand=0; iVertexCand<nVertexCand_; ++iVertexCand )
  {
    vertexMass_ = vertexPt_ = -1e9;
    vertexL3D_ = vertexL2D_ = vertexLxy_ = -1e9;
    deltaR_ = deltaPt_ = 1e-9;
    legId_ = -1;
    muonCharge_ = trackCharge_ = 0;
    muon_ = track_ = math::XYZTLorentzVector();

    const reco::VertexCompositeCandidate& vertexCand = vertexCandHandle->at(iVertexCand);
    if ( !(*vertexCut_)(vertexCand) ) continue;

    vertexMass_ = vertexCand.mass();
    vertexPt_ = vertexCand.pt();

    math::XYZPoint vertex = vertexCand.vertex();
    vertexL3D_ = vertex.R();
    vertexL2D_ = vertex.Rho();
    vertexLxy_ = (vertex.x()*vertexCand.px() + vertex.y()*vertexCand.py())/vertexCand.p();

    const reco::Candidate* p1 = vertexCand.daughter(0);
    const reco::Candidate* p2 = vertexCand.daughter(1);

    // Count duplicated track
    bool isOverlapChecked = false;
    for ( int j=0, n=overlappedCandIds.size(); j<n; ++j )
    {
      if ( iVertexCand != overlappedCandIds[j] ) continue;

      isOverlapChecked = true;
      break;
    }

    if ( !isOverlapChecked )
    {
      int nOverlap = 0;
      for ( int jVertexCand=iVertexCand+1; jVertexCand<nVertexCand_; ++jVertexCand )
      {
        const reco::VertexCompositeCandidate& vertexCand2 = vertexCandHandle->at(jVertexCand);
        const reco::Candidate* pp1 = vertexCand2.daughter(0);
        const reco::Candidate* pp2 = vertexCand2.daughter(1);

        if ( p1->get<reco::TrackRef>() == pp1->get<reco::TrackRef>() or
             p2->get<reco::TrackRef>() == pp2->get<reco::TrackRef>() or
             p1->get<reco::TrackRef>() == pp2->get<reco::TrackRef>() or
             p2->get<reco::TrackRef>() == pp1->get<reco::TrackRef>() )
        {
          // Detected overlap
          ++nOverlap;
        }
      }
      hNCandOverlap_->Fill(nOverlap);
      overlappedCandIds.push_back(iVertexCand);
    }

    // Find fake muon
    const reco::Muon* muon1 = findMatchedMuonByTrackRef(*p1, muonHandle);
    const reco::Muon* muon2 = findMatchedMuonByTrackRef(*p2, muonHandle);

    const reco::Muon* matchedMuon = muon1 ? muon1 : muon2;
    const reco::Candidate* matchedTrack = muon1 ? p1 : p2;
    // Set matched muon. 
    // legId = 0 for un-matched case.
    // There can be very rare case of double fake muons. 
    // Only muon1 will be stored in this case but we can separate this case with requiring legId != 3
    legId_ = 0;
    if ( matchedMuon )
    {
      if ( muon1 ) { legId_ |= 1; ++nFakeMuon; }
      if ( muon2 ) { legId_ |= 2; ++nFakeMuon; }
      
      muonCharge_ = matchedMuon->charge();
      trackCharge_ = matchedTrack->charge();
      muon_ = matchedMuon->p4();
      track_ = matchedTrack->p4();
      deltaR_ = deltaR(*matchedTrack, *matchedMuon);
      deltaPt_ = matchedTrack->pt() - matchedMuon->pt();

      for ( int i=0, n=muonTypes_.size(); i<n; ++i )
      {
        muonTypeResults_[i] = (*muonTypes_[i])(*matchedMuon);
      }
    }

    tree_->Fill();
  }

  if ( nFakeMuon != 0 ) hEvent_->Fill(3);
}

const reco::Muon* FakeMuonAnalyzer::findMatchedMuonByTrackRef(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle)
{
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const reco::Muon& muonCand = muonHandle->at(i);

    if ( p.get<reco::TrackRef>() == muonCand.track() ) return &muonCand;
  }

  return 0;
}

const reco::Muon* FakeMuonAnalyzer::findMatchedMuonByDR(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle)
{
  const reco::Muon* matchedMuon = 0;
  double matchedDR = maxDR_;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const reco::Muon& muonCand = muonHandle->at(i);

    const double dR = deltaR(p, muonCand);
    if ( dR > matchedDR ) continue;

    matchedDR = dR;
    matchedMuon = &muonCand;
  }

  return matchedMuon;
}

const reco::Muon* FakeMuonAnalyzer::findMatchedMuonByDRDPt(const reco::Candidate& p, edm::Handle<edm::View<reco::Muon> >& muonHandle)
{
  const reco::Muon* matchedMuon = 0;
  double matchedDR = maxDR_;
  double matchedDPt = maxDPt_;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const reco::Muon& muonCand = muonHandle->at(i);

    const double dR = deltaR(p, muonCand);
    if ( dR > matchedDR ) continue;
    const double dPt = abs(p.pt() - muonCand.pt());
    if ( dPt > matchedDPt ) continue;

    matchedDR = dR;
    matchedDPt = dPt;
    matchedMuon = &muonCand;
  }

  return matchedMuon;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FakeMuonAnalyzer);
