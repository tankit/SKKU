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

  Int_t runNumber, eventNumber, nMuon, nSelMuon;
  Double_t muPt, muP, muEta, muPhi;
  Bool_t trkMu, trkMuArb, rpcMu, rpcMuLoose, rpcMuMedium, rpcMuTight;
  Bool_t staMu, glbMu, glbMuPromptT, glbMuLoose, glbMuMedium, glbMuTight;

  TH1F* hNMuon_;
  TH1F* hNRPCMuon_;
  TH1F* hNGlbPromptT_;
  TH1F* hNRPCMuMedium_;
  TH1F* hNTrkArbitrated_;
  TH1F* hNRPCMuLoose_;
  TH1F* hNRPCMuTight_;
  TH1F* hNGlbMuTight_;

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
  tree_->Branch("muPt",         &muPt,         "muPt/D");
  tree_->Branch("muP",          &muP,          "muP/D");
  tree_->Branch("muEta",        &muEta,        "muEta/D");
  tree_->Branch("muPhi",        &muPhi,        "muPhi/D");
  tree_->Branch("trkMu",        &trkMu,        "trkMu/B");
  tree_->Branch("trkMuArb",     &trkMuArb,     "trkMuArb/B");
  tree_->Branch("rpcMu",        &rpcMu,        "rpcMu/B");
  tree_->Branch("rpcMuLoose",   &rpcMuLoose,   "rpcMuLoose/B");
  tree_->Branch("rpcMuMedium",  &rpcMuMedium,  "rpcMuMedium/B");
  tree_->Branch("rpcMuTight",   &rpcMuTight,   "rpcMuTight/B");
  tree_->Branch("staMu",        &staMu,        "staMu/B");
  tree_->Branch("glbMu",        &glbMu,        "glbMu/B");
  tree_->Branch("glbMuLoose",   &glbMuLoose,   "glbMuLoose/B");
  tree_->Branch("glbMuPromptT", &glbMuPromptT, "glbMuPromptT/B");
  tree_->Branch("glbMuMedium",  &glbMuMedium,  "glbMuMedium/B");
  tree_->Branch("glbMuTight",   &glbMuTight,   "glbMuTight/B");

  hNMuon_       = fs->make<TH1F>("hNMuon"      , "Number of muons;Number of muons", 10, 0, 10);
  hNRPCMuon_    = fs->make<TH1F>("hNRPCMuon"   , "Number of RPC muons;Number of muons", 10, 0, 10);
  hNGlbPromptT_ = fs->make<TH1F>("hNGlbPromptT", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNRPCMuMedium_ = fs->make<TH1F>("hNRPCMuMedium", "Number of RPCMuMedium muons;Number of muons", 10, 0, 10);
  hNTrkArbitrated_ = fs->make<TH1F>("hNTrkArbitrated", "Number of TrkMuArbitrated muons;Number of muons", 10, 0, 10);
  hNRPCMuLoose_  = fs->make<TH1F>("hNRPCMuLoose"   , "Number of RPC muons;Number of muons", 10, 0, 10);
  hNRPCMuTight_  = fs->make<TH1F>("hNRPCMuTight"   , "Number of RPC muons;Number of muons", 10, 0, 10);
  hNGlbMuTight_  = fs->make<TH1F>("hNGlbMuTight"   , "Number of RPC muons;Number of muons", 10, 0, 10);

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
  
  nMuon = muonHandle->size(); nSelMuon = 0;
  int nGlbPromptT = 0, nTrkArbitrated = 0, nGlbMuTight = 0;
  int nRPCMuon = 0, nRPCMuMedium = 0, nRPCMuLoose = 0, nRPCMuTight = 0;
  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {

    glbMu = staMu = trkMu = rpcMu = rpcMuLoose = rpcMuMedium = rpcMuTight = trkMuArb = glbMuLoose = glbMuPromptT = glbMuMedium = glbMuTight = false; 
    if ( muon->pt() < minPtTrk_ ) continue; 
    const double abseta = abs(muon->eta());
    if ( abseta > maxEtaTrk_ ) continue;

    muPt = muon->pt();
    muP  = muon->p();
    muEta = muon->eta();
    muPhi = muon->phi(); 
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
    //--Note for Z MC: 100% efficient with muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 because of "GlobalMuLoose = GlobalMuonPromptTight" and "GlobalMuMedium = GlobalMuTight"

    ++nSelMuon;
    if ( idFlags[1] ) glbMu = true;
    if ( idFlags[2] ) staMu = true;
    if ( idFlags[3] ) trkMu = true;

    if ( idFlags[4] )
    {
      ++nRPCMuon; rpcMu = true;

      if ( idFlags[5] ) { ++nRPCMuLoose;  rpcMuLoose  = true; }
      if ( idFlags[6] ) { ++nRPCMuMedium; rpcMuMedium = true; }
      if ( idFlags[7] ) { ++nRPCMuTight;  rpcMuTight  = true; }
    }

    if ( idFlags[8] )  { ++nTrkArbitrated; trkMuArb = true; }
    if ( idFlags[9] ) glbMuLoose = true;
    if ( idFlags[10] ) { ++nGlbPromptT;    glbMuPromptT = true; }
    if ( idFlags[11] ) glbMuMedium = true;
    if ( idFlags[12] ) { ++nGlbMuTight; glbMuTight = true; }

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
  hNTrkArbitrated_->Fill(nTrkArbitrated);
  hNRPCMuLoose_->Fill(nRPCMuLoose);
  hNRPCMuTight_->Fill(nRPCMuTight);
  hNGlbMuTight_->Fill(nGlbMuTight);

  tree_->Fill();

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonAnalyzer);
