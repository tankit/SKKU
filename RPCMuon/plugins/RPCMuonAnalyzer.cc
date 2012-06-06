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
  std::vector<bool> trkMuArb;
  std::vector<bool> rpcMu;
  std::vector<bool> rpcMuLoose;
  std::vector<bool> rpcMuMedium;
  std::vector<bool> rpcMuTight;
  std::vector<bool> staMu;
  std::vector<bool> glbMu;
  std::vector<bool> glbMuLoose;
  std::vector<bool> glbMuPromptT;
  std::vector<bool> glbMuMedium;
  std::vector<bool> glbMuTight;

  Int_t runNumber, eventNumber, nMuon, nSelMuon;
  Int_t nGlbMuon, nStaMuon, nTrkMuon;
  Int_t nGlbPromptT, nTrkMuArb, nGlbMuLoose, nGlbMuMedium, nGlbMuTight;
  Int_t nRPCMuon, nRPCMuMedium, nRPCMuLoose, nRPCMuTight;

  TH1F* hNMuon_;
  TH1F* hNRPCMuon_;
  TH1F* hNGlbPromptT_;
  TH1F* hNRPCMuMedium_;
  TH1F* hNTrkMuArb_;
  TH1F* hNRPCMuLoose_;
  TH1F* hNRPCMuTight_;
  TH1F* hNGlbMuTight_;
  TH1F* hNTrkMuon_;

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
  tree_->Branch("nRPCMuLoose",  &nRPCMuLoose,  "nRPCMuLoose/I");
  tree_->Branch("nRPCMuMedium", &nRPCMuMedium, "nRPCMuMedium/I");
  tree_->Branch("nRPCMuTight",  &nRPCMuTight,  "nRPCMuTight/I");
  tree_->Branch("nStaMuon",     &nStaMuon,     "nStaMuon/I");
  tree_->Branch("nTrkMuon",     &nTrkMuon,     "nTrkMuon/I");
  tree_->Branch("nTrkMuArb",    &nTrkMuArb,    "nTrkMuArb/I");
  tree_->Branch("nGlbMuon",     &nGlbMuon,     "nGlbMuon/I");
  tree_->Branch("nGlbMuLoose",  &nGlbMuLoose,  "nGlbMuLoose/I");
  tree_->Branch("nGlbMuMedium", &nGlbMuMedium, "nGlbMuMedium/I");
  tree_->Branch("nGlbMuTight",  &nGlbMuTight,  "nGlbMuTight/I");
  tree_->Branch("nGlbPromptT",  &nGlbPromptT,  "nGlbPromptT/I");
  tree_->Branch("muPt",         &muPt);
  tree_->Branch("muP",          &muP);
  tree_->Branch("muEta",        &muEta);
  tree_->Branch("muPhi",        &muPhi);
  tree_->Branch("trkMu",        &trkMu);
  tree_->Branch("trkMuArb",     &trkMuArb);
  tree_->Branch("rpcMu",        &rpcMu);
  tree_->Branch("rpcMuLoose",   &rpcMuLoose);
  tree_->Branch("rpcMuMedium",  &rpcMuMedium);
  tree_->Branch("rpcMuTight",   &rpcMuTight);
  tree_->Branch("staMu",        &staMu);
  tree_->Branch("glbMu",        &glbMu);
  tree_->Branch("glbMuLoose",   &glbMuLoose);
  tree_->Branch("glbMuPromptT", &glbMuPromptT);
  tree_->Branch("glbMuMedium",  &glbMuMedium);
  tree_->Branch("glbMuTight",   &glbMuTight);

  hNMuon_       = fs->make<TH1F>("hNMuon"      , "Number of muons;Number of muons", 10, 0, 10);
  hNRPCMuon_    = fs->make<TH1F>("hNRPCMuon"   , "Number of RPC muons;Number of muons", 10, 0, 10);
  hNGlbPromptT_ = fs->make<TH1F>("hNGlbPromptT", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNRPCMuMedium_ = fs->make<TH1F>("hNRPCMuMedium", "Number of RPCMuMedium muons;Number of muons", 10, 0, 10);
  hNTrkMuArb_    = fs->make<TH1F>("hNTrkMuArb"  , "Number of TrkMuArbitrated muons;Number of muons", 10, 0, 10);
  hNRPCMuLoose_  = fs->make<TH1F>("hNRPCMuLoose", "Number of RPCMuLoose;Number of muons", 10, 0, 10);
  hNRPCMuTight_  = fs->make<TH1F>("hNRPCMuTight", "Number of RPCMuTight;Number of muons", 10, 0, 10);
  hNGlbMuTight_  = fs->make<TH1F>("hNGlbMuTight", "Number of GlbMuTight;Number of muons", 10, 0, 10);
  hNTrkMuon_     = fs->make<TH1F>("hNTrkMuon"   , "Number of Tracker muons;Number of muons", 10, 0, 10);

  const char* idNames[] = {
    "All", "AllGlbMu", "AllStaMu", "AllTrkMu", "AllRPCMu", "RPCMuLoose", "RPCMuMedium", "RPCMuTight", "TrkMuArbitrated", "GlbMuLoose", "GlbPromptTight", "GlbMuMedium", "GlbMuTight"
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
  trkMuArb.clear();
  rpcMu.clear();
  rpcMuLoose.clear();
  rpcMuMedium.clear();
  rpcMuTight.clear();
  staMu.clear();
  glbMu.clear();
  glbMuLoose.clear();
  glbMuPromptT.clear();
  glbMuMedium.clear();
  glbMuTight.clear();

  nMuon = muonHandle->size(); nSelMuon = 0;
  nGlbMuon = 0, nStaMuon = 0, nTrkMuon = 0;
  nGlbPromptT = 0, nTrkMuArb = 0, nGlbMuLoose = 0, nGlbMuMedium = 0, nGlbMuTight = 0;
  nRPCMuon = 0, nRPCMuMedium = 0, nRPCMuLoose = 0, nRPCMuTight = 0;
  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {

    if ( muon->pt() < minPtTrk_ ) continue; 
    const double abseta = abs(muon->eta());
    if ( abseta > maxEtaTrk_ ) continue;

    std::cout << " * Muon Pt = " << muon->pt() << ", P = " << muon->p() << ", Eta = " << muon->eta() << ", Phi = " << muon->phi() << std::endl;

    const bool idFlags[] = {
      true,
      muon->isGlobalMuon(), muon->isStandAloneMuon(), muon->isTrackerMuon(),
      muon->isRPCMuon(),
      muon::isGoodMuon(*muon, muon::RPCMuLoose),
      muon::isGoodMuon(*muon, muon::RPCMuLoose) && muon->numberOfMatchedStations(reco::Muon::RPCHitAndTrackArbitration)>1,
      muon->isRPCMuon() && muon::isGoodMuon(*muon, muon::RPCMu, 3, 20, 4, 1e9, 1e9, 1e9, 1e9, reco::Muon::RPCHitAndTrackArbitration, false, false),
      muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated),
      muon->isGlobalMuon() && muon->globalTrack()->normalizedChi2()<10., 
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight),
      muon->isGlobalMuon() && muon->globalTrack()->normalizedChi2()<10. && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1,
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight) && muon->numberOfMatchedStations(reco::Muon::SegmentAndTrackArbitration)>1,

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

      if ( idFlags[5] ) ++nRPCMuLoose;
      if ( idFlags[6] ) ++nRPCMuMedium;
      if ( idFlags[7] ) ++nRPCMuTight;
    }

    if ( idFlags[8] ) ++nTrkMuArb;
    if ( idFlags[9] ) ++nGlbMuLoose;
    if ( idFlags[10] ) ++nGlbPromptT;
    if ( idFlags[11] ) ++nGlbMuMedium;
    if ( idFlags[12] ) ++nGlbMuTight;

    muPt.push_back(muon->pt());
    muP.push_back(muon->p());
    muEta.push_back(muon->eta());
    muPhi.push_back(muon->phi());

    glbMu.push_back(idFlags[1]);
    staMu.push_back(idFlags[2]);
    trkMu.push_back(idFlags[3]);
    rpcMu.push_back(idFlags[4]);
    rpcMuLoose.push_back(idFlags[5]);
    rpcMuMedium.push_back(idFlags[6]);
    rpcMuTight.push_back(idFlags[7]);
    trkMuArb.push_back(idFlags[8]);
    glbMuLoose.push_back(idFlags[9]);
    glbMuPromptT.push_back(idFlags[10]);
    glbMuMedium.push_back(idFlags[11]);
    glbMuTight.push_back(idFlags[12]);

    std::cout << " + idFlags [RPCMu, Loose, Medium, Tight] = " << idFlags[4] << " " << idFlags[5] << " " << idFlags[6] << " " << idFlags[7] << std::endl;

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
  hNRPCMuMedium_->Fill(nRPCMuMedium);
  hNGlbPromptT_->Fill(nGlbPromptT);
  hNTrkMuArb_->Fill(nTrkMuArb);
  hNRPCMuLoose_->Fill(nRPCMuLoose);
  hNRPCMuTight_->Fill(nRPCMuTight);
  hNGlbMuTight_->Fill(nGlbMuTight);
  hNTrkMuon_->Fill(nTrkMuon);

  tree_->Fill();

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonAnalyzer);
