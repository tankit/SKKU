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

  TTree* tree_;

  TH1F* hNMuon_;
  TH1F* hNRPCMuon_;
  TH1F* hNGlbPromptT_;
  TH1F* hNRPCMuBasic_;

  TH2F* hIdCorrelation_;
  TH2F* hBarrel_IdCorrelation_;
  TH2F* hEndcap_IdCorrelation_;
};

RPCMuonAnalyzer::RPCMuonAnalyzer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getUntrackedParameter<edm::InputTag>("muon");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");

  hNMuon_       = fs->make<TH1F>("hNMuon"      , "Number of muons;Number of muons", 10, 0, 10);
  hNRPCMuon_    = fs->make<TH1F>("hNRPCMuon"   , "Number of RPC muons;Number of muons", 10, 0, 10);
  hNGlbPromptT_ = fs->make<TH1F>("hNGlbPromptT", "Number of GlobalMuPromptTight muons;Number of muons", 10, 0, 10);
  hNRPCMuBasic_ = fs->make<TH1F>("hNRPCMuBasic", "Number of RPCMuBasic muons;Number of muons", 10, 0, 10);

  const char* idNames[] = {
    "All", "Glb", "Sta", "Trk", "RPC", "GlbPromptT", "RPCBasic", "RPCLoose", "RPCVeryLoose"
  };
  const int nId = sizeof(idNames)/sizeof(const char*);
  hIdCorrelation_ = fs->make<TH2F>("hIdCorrelation", "ID correlation", nId, 0, nId, nId, 0, nId);
  hBarrel_IdCorrelation_ = fs->make<TH2F>("hIdCorrelationBarrel", "ID correlation (Barrel)", nId, 0, nId, nId, 0, nId);
  hEndcap_IdCorrelation_ = fs->make<TH2F>("hIdCorrelationEndcap", "ID correlation (Endcap)", nId, 0, nId, nId, 0, nId);
  for ( int i=0; i<nId; ++i )
  {
    hIdCorrelation_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hIdCorrelation_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
    hBarrel_IdCorrelation_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hBarrel_IdCorrelation_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
    hEndcap_IdCorrelation_->GetXaxis()->SetBinLabel(i+1, idNames[i]);
    hEndcap_IdCorrelation_->GetYaxis()->SetBinLabel(i+1, idNames[i]);
  }
  hIdCorrelation_->SetOption("COLZ");
  hBarrel_IdCorrelation_->SetOption("COLZ");
  hEndcap_IdCorrelation_->SetOption("COLZ");
}

void RPCMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  
  int nMuon = muonHandle->size();
  int nGlbPromptT = 0;
  int nRPCMuon = 0, nRPCMuBasic = 0;
  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {
    if ( muon->pt() < 5 ) continue;
    const double abseta = abs(muon->eta());
    if ( abseta > 2.1 ) continue;

    const bool idFlags[] = {
      true,
      muon->isGlobalMuon(), muon->isStandAloneMuon(), muon->isTrackerMuon(),
      muon->isRPCMuon(), 
      muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight),
      muon::isGoodMuon(*muon, muon::RPCMuBasic),
      muon::isGoodMuon(*muon, muon::RPCMuLoose),
      muon::isGoodMuon(*muon, muon::RPCMuVeryLoose),
    };

    if ( idFlags[4] )
    {
      ++nRPCMuon;

      if ( idFlags[6] ) ++nRPCMuBasic;
    }

    if ( idFlags[5] ) ++nGlbPromptT;

    // Fill correlation matrix
    for ( int i=0, n=sizeof(idFlags)/sizeof(const bool); i<n; ++i )
    {
      for ( int j=i; j<n; ++j )
      {
        if ( idFlags[i] and idFlags[j] )
        {
          hIdCorrelation_->Fill(i,j);
          if ( abseta < 1.6 ) hBarrel_IdCorrelation_->Fill(i,j);
          else hEndcap_IdCorrelation_->Fill(i,j);
        }
      }
    }
  }

  hNMuon_->Fill(nMuon);
  hNRPCMuon_->Fill(nRPCMuon);
  hNRPCMuBasic_->Fill(nRPCMuBasic);
  hNGlbPromptT_->Fill(nGlbPromptT);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonAnalyzer);
