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
  TH1F* hNRPCMuBasic_;

};

RPCMuonAnalyzer::RPCMuonAnalyzer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getUntrackedParameter<edm::InputTag>("muon");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
  hNMuon_ = fs->make<TH1F>("hNMuon", "Number of muons;Number of muons", 10, 0, 10);
  hNRPCMuon_ = fs->make<TH1F>("hNRPCMuon", "NUmber of RPC muons;Number of muons", 10, 0, 10);
  hNRPCMuBasic_ = fs->make<TH1F>("hNRPCMuBasic", "Number of RPCMuBasic muons;Number of muons", 10, 0, 10);
}

void RPCMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  
  int nMuon = muonHandle->size();
  int nRPCMuon = 0, nRPCMuBasic = 0;
  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {
    if ( !muon->isRPCMuon() ) continue;
    ++nRPCMuon;

    if ( !muon::isGoodMuon(*muon, muon::RPCMuBasic) ) continue;
    ++nRPCMuBasic;
  }

  hNMuon_->Fill(nMuon);
  hNRPCMuon_->Fill(nRPCMuon);
  hNRPCMuBasic_->Fill(nRPCMuBasic);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonAnalyzer);
