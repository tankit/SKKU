#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "DataFormats/TrackReco/interface/Track.h"
//#include "Geometry/Records/interface/MuonGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

class RPCMuonAnalyzer : public edm::EDAnalyzer 
{
public:
  RPCMuonAnalyzer(const edm::ParameterSet& pset);
  ~RPCMuonAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  void beginJob() {};
  void endJob() {};

private:
  edm::InputTag muonLabel_;
  double minPtTrk_;
  double maxEtaTrk_;

  TTree* tree_;

  std::vector<float> muPt;
  std::vector<float> muP;
  std::vector<float> muEta;
  std::vector<float> muPhi;
  std::vector<bool> trkMu;
  std::vector<bool> trkMuTight;
  std::vector<bool> trkMuTight2;
  std::vector<bool> rpcMu;
  std::vector<bool> rpcMuTight;
  std::vector<bool> staMu;
  std::vector<bool> glbMu;
  std::vector<bool> glbMuTight;
  std::vector<bool> glbMuTight2;
  std::vector<bool> glbMuTighter;
  std::vector<bool> glbMuTighter2;
  std::vector<int> nMatchesSegNoArb;
  std::vector<int> nStationSegNoArb;
  std::vector<int> nMatchesSegTrkArb;
  std::vector<int> nStationSegTrkArb;
  std::vector<int> nMatchesRPCTrkArb;
  std::vector<int> nStationRPCTrkArb;
  std::vector<int> nLayerRPCTrkArb;

  Int_t runNumber, eventNumber, nMuon, nSelMuon;
  Int_t nGlbMuon, nStaMuon, nTrkMuon;
  Int_t nRPCMuon, nRPCMuTight;
  Int_t nTrkMuTight, nTrkMuTight2;
  Int_t nGlbMuTight, nGlbMuTight2, nGlbMuTighter, nGlbMuTighter2;

  TH1F* hNMuon_;
  TH1F* hNRPCMuon_;
  TH1F* hNRPCMuTight_;
  TH1F* hNTrkMuTight_;
  TH1F* hNTrkMuTight2_;
  TH1F* hNGlbMuTight_;
  TH1F* hNGlbMuTight2_;
  TH1F* hNGlbMuTighter_;
  TH1F* hNGlbMuTighter2_;

  TH2F* hIdCorrelation_;
  TH2F* hIdCorrelationB_;
  TH2F* hIdCorrelationO_;
  TH2F* hIdCorrelationE_;
};

RPCMuonAnalyzer::RPCMuonAnalyzer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getUntrackedParameter<edm::InputTag>("muon");
  minPtTrk_  = pset.getUntrackedParameter<double>("minPtTrk");
  maxEtaTrk_ = pset.getUntrackedParameter<double>("maxEtaTrk");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("runNumber",    &runNumber,    "runNumber/I");
  tree_->Branch("eventNumber",  &eventNumber,  "eventNumber/I");
  tree_->Branch("nMuon",        &nMuon,        "nMuon/I");
  tree_->Branch("nSelMuon",     &nSelMuon,     "nSelMuon/I");
  tree_->Branch("nRPCMuon",     &nRPCMuon,     "nRPCMuon/I");
  tree_->Branch("nRPCMuTight",  &nRPCMuTight,  "nRPCMuTight/I");
  tree_->Branch("nStaMuon",     &nStaMuon,     "nStaMuon/I");
  tree_->Branch("nTrkMuon",     &nTrkMuon,     "nTrkMuon/I");
  tree_->Branch("nGlbMuon",     &nGlbMuon,     "nGlbMuon/I");
  tree_->Branch("nTrkMuTight",  &nTrkMuTight,  "nTrkMuTight/I");
  tree_->Branch("nTrkMuTight2", &nTrkMuTight2, "nTrkMuTight2/I");  
  tree_->Branch("nGlbMuTight",     &nGlbMuTight,     "nGlbMuTight/I");
  tree_->Branch("nGlbMuTight2",    &nGlbMuTight2,    "nGlbMuTight2/I");
  tree_->Branch("nGlbMuTighter",   &nGlbMuTighter,   "nGlbMuTighter/I");
  tree_->Branch("nGlbMuTighter2",  &nGlbMuTighter2,  "nGlbMuTighter2/I");
  tree_->Branch("muPt",         &muPt);
  tree_->Branch("muP",          &muP);
  tree_->Branch("muEta",        &muEta);
  tree_->Branch("muPhi",        &muPhi);
  tree_->Branch("trkMu",        &trkMu);
  tree_->Branch("trkMuTight",   &trkMuTight);
  tree_->Branch("trkMuTight2",  &trkMuTight2);
  tree_->Branch("rpcMu",        &rpcMu);
  tree_->Branch("rpcMuTight",   &rpcMuTight);
  tree_->Branch("staMu",        &staMu);
  tree_->Branch("glbMu",        &glbMu);
  tree_->Branch("glbMuTight",    &glbMuTight);
  tree_->Branch("glbMuTight2",   &glbMuTight2);
  tree_->Branch("glbMuTighter",  &glbMuTighter);
  tree_->Branch("glbMuTighter2", &glbMuTighter2);
  tree_->Branch("nMatchesSegNoArb",  &nMatchesSegNoArb);
  tree_->Branch("nStationSegNoArb",  &nStationSegNoArb);
  tree_->Branch("nMatchesSegTrkArb", &nMatchesSegTrkArb);
  tree_->Branch("nStationSegTrkArb", &nStationSegTrkArb);
  tree_->Branch("nMatchesRPCTrkArb", &nMatchesRPCTrkArb);
  tree_->Branch("nStationRPCTrkArb", &nStationRPCTrkArb);
  tree_->Branch("nLayerRPCTrkArb",   &nLayerRPCTrkArb);

  hNMuon_          = fs->make<TH1F>("hNMuon", "Number of muons;Number of muons", 10, 0, 10);
  hNRPCMuon_       = fs->make<TH1F>("hNRPCMuon", "Number of RPC muons;Number of muons", 10, 0, 10);
  hNRPCMuTight_    = fs->make<TH1F>("hNRPCMuTight", "Number of RPCMuTight;Number of muons", 10, 0, 10);
  hNTrkMuTight_    = fs->make<TH1F>("hNTrkMuTight", "Number of TrkMuTight muons;Number of muons", 10, 0, 10);
  hNTrkMuTight2_   = fs->make<TH1F>("hNTrkMuTight2", "Number of TrkMuTight muons;Number of muons", 10, 0, 10);
  hNGlbMuTight_    = fs->make<TH1F>("hNGlbMuTight", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNGlbMuTight2_   = fs->make<TH1F>("hNGlbMuTight2", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNGlbMuTighter_  = fs->make<TH1F>("hNGlbMuTighter", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNGlbMuTighter2_ = fs->make<TH1F>("hNGlbMuTighter2", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);

  const char* idNames[] = {
    "All", "AllGlbMu", "AllStaMu", "AllTrkMu", "AllRPCMu", "RPCMuTight", "TMOneStationTight", "TMOneStationTight+", "GlbPromptTight", "GlbPromptTight+", "GlbPromptTighter", "GlbPromptTighter+"
  };
  const int nId = sizeof(idNames)/sizeof(const char*);
  hIdCorrelation_ = fs->make<TH2F>("hIdCorrelation", "ID correlation", nId, 0, nId, nId, 0, nId);
  hIdCorrelationB_ = fs->make<TH2F>("hIdCorrelationBarrel", "ID correlation (Barrel)", nId, 0, nId, nId, 0, nId);
  hIdCorrelationO_ = fs->make<TH2F>("hIdCorrelationOverlap", "ID correlation (Overlap)", nId, 0, nId, nId, 0, nId);
  hIdCorrelationE_ = fs->make<TH2F>("hIdCorrelationEndcap", "ID correlation (Endcap)", nId, 0, nId, nId, 0, nId);
  for ( int i=0; i<nId; ++i )
  {
    hIdCorrelation_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelation_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationB_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationB_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationO_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationO_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationE_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelationE_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
  }
  hIdCorrelation_->SetOption("COLZ");
  hIdCorrelationB_->SetOption("COLZ");
  hIdCorrelationO_->SetOption("COLZ");
  hIdCorrelationE_->SetOption("COLZ");
}

void RPCMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{

  // select the event
  runNumber   = event.id().run();
  eventNumber = event.id().event();

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  
  muPt.clear();
  muP.clear();
  muEta.clear();
  muPhi.clear();

  trkMu.clear();
  trkMuTight.clear();
  trkMuTight2.clear();
  rpcMu.clear();
  rpcMuTight.clear();
  staMu.clear();
  glbMu.clear();
  glbMuTight.clear();
  glbMuTight2.clear();
  glbMuTighter.clear();
  glbMuTighter2.clear();

  nMatchesSegNoArb.clear();
  nStationSegNoArb.clear();
  nMatchesSegTrkArb.clear();
  nStationSegTrkArb.clear();
  nMatchesRPCTrkArb.clear();
  nStationRPCTrkArb.clear();
  nLayerRPCTrkArb.clear();

  nMuon = muonHandle->size(); nSelMuon = 0;
  nGlbMuon = 0, nStaMuon = 0, nTrkMuon = 0;
  nRPCMuon = 0, nRPCMuTight = 0;
  nTrkMuTight = 0, nTrkMuTight2 = 0;
  nGlbMuTight = 0, nGlbMuTight2 = 0, nGlbMuTighter = 0, nGlbMuTighter2 = 0;
  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {

    if ( muon->pt() < minPtTrk_ ) continue; 
    const double abseta = abs(muon->eta());
    if ( abseta > maxEtaTrk_ ) continue;

    //std::cout << " * Muon Pt = " << muon->pt() << ", P = " << muon->p() << ", Eta = " << muon->eta() << ", Phi = " << muon->phi() << std::endl;

    const bool idFlags[] = {
      true,
      muon->isGlobalMuon(), muon->isStandAloneMuon(), muon->isTrackerMuon(),
      muon->isRPCMuon(),
      muon::isGoodMuon(*muon, muon::RPCMuLoose) && muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration)>1 && muon->numberOfMatchedLayers(reco::Muon::RPCHitAndTrackArbitration)>2,
      muon::isGoodMuon(*muon, muon::TMOneStationTight),
      muon::isGoodMuon(*muon, muon::TMOneStationTight) || (muon::isGoodMuon(*muon, muon::RPCMuLoose) && muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration)>1 && muon->numberOfMatchedLayers(reco::Muon::RPCHitAndTrackArbitration)>2),
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight),
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight) || (muon::isGoodMuon(*muon, muon::RPCMuLoose) && muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration)>1 && muon->numberOfMatchedLayers(reco::Muon::RPCHitAndTrackArbitration)>2),
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight) && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1,
      (muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight) && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1) || (muon::isGoodMuon(*muon, muon::RPCMuLoose) && muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration)>1 && muon->numberOfMatchedLayers(reco::Muon::RPCHitAndTrackArbitration)>2)
    };

    //--GlobalMuLoose
    //muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10.
    //--GlobalMuonPromptTight
    //muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0
    //--GlobalMuMedium
    //muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1
    //--GlobalMuTight
    //GlobalMuMedium && muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0
    // or
    //GlobalMuonPromptTight && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1
    //
    //--Note for Z MC: how efficient with muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0? e.g, "GlobalMuLoose vs. GlobalMuonPromptTight" and "GlobalMuMedium vs. GlobalMuTight"

    ++nSelMuon;
    if ( idFlags[1] ) ++nGlbMuon;
    if ( idFlags[2] ) ++nStaMuon;
    if ( idFlags[3] ) ++nTrkMuon;

    if ( idFlags[4] )
    {
      ++nRPCMuon;

      if ( idFlags[5] ) ++nRPCMuTight;
      if ( idFlags[6] ) ++nTrkMuTight;
      if ( idFlags[7] ) ++nTrkMuTight2;
    }

    if ( idFlags[8] ) ++nGlbMuTight;
    if ( idFlags[9] ) ++nGlbMuTight2;
    if ( idFlags[10] ) ++nGlbMuTighter;
    if ( idFlags[11] ) ++nGlbMuTighter2;

    muPt.push_back(muon->pt());
    muP.push_back(muon->p());
    muEta.push_back(muon->eta());
    muPhi.push_back(muon->phi());

    glbMu.push_back(idFlags[1]);
    staMu.push_back(idFlags[2]);
    trkMu.push_back(idFlags[3]);
    rpcMu.push_back(idFlags[4]);
    rpcMuTight.push_back(idFlags[5]);
    trkMuTight.push_back(idFlags[6]);
    trkMuTight2.push_back(idFlags[7]);
    glbMuTight.push_back(idFlags[8]);
    glbMuTight2.push_back(idFlags[9]);
    glbMuTighter.push_back(idFlags[10]);
    glbMuTighter2.push_back(idFlags[11]);

    nMatchesSegNoArb.push_back(muon->numberOfMatches(reco::Muon::NoArbitration));
    nStationSegNoArb.push_back(muon->numberOfMatchedStations(reco::Muon::NoArbitration));
    nMatchesSegTrkArb.push_back(muon->numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
    nStationSegTrkArb.push_back(muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration));
    nMatchesRPCTrkArb.push_back(muon->numberOfMatches(reco::Muon::RPCHitAndTrackArbitration));
    nStationRPCTrkArb.push_back(muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration));
    nLayerRPCTrkArb.push_back(muon->numberOfMatchedLayers(reco::Muon::RPCHitAndTrackArbitration));

    //std::cout << " + idFlags [RPCMu, Loose, Medium, Tight] = " << idFlags[4] << " " << idFlags[5] << " " << idFlags[6] << " " << idFlags[7] << std::endl;

    // Fill correlation matrix
    for ( int i=0, n=sizeof(idFlags)/sizeof(const bool); i<n; ++i )
    {
      for ( int j=i; j<n; ++j )
      {
        if ( idFlags[i] and idFlags[j] )
        {
          hIdCorrelation_->Fill(i,j);
          if ( abseta < 0.8 ) hIdCorrelationB_->Fill(i,j);
          else if ( abseta < 1.2 ) hIdCorrelationO_->Fill(i,j);
          else hIdCorrelationE_->Fill(i,j);
        }
      }
    }
  }

  hNMuon_->Fill(nMuon);
  hNRPCMuon_->Fill(nRPCMuon);
  hNRPCMuTight_->Fill(nRPCMuTight);
  hNTrkMuTight_->Fill(nTrkMuTight);
  hNTrkMuTight2_->Fill(nTrkMuTight2);
  hNGlbMuTight_->Fill(nGlbMuTight);
  hNGlbMuTight2_->Fill(nGlbMuTight2);
  hNGlbMuTighter_->Fill(nGlbMuTighter);
  hNGlbMuTighter2_->Fill(nGlbMuTighter2);

  tree_->Fill();

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonAnalyzer);
