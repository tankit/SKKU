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
  edm::InputTag muonLabel_;
  edm::InputTag kshortLabel_;
  double maxDR_;
  std::vector<StringCutObjectSelector<reco::Muon, true>*> idCuts_;

  TH1F* hEvent_;
  TTree* tree_;
  int run_, lumi_, event_;
  int nMuon_, nKshort_;
  std::vector<double>* masses_;
  std::vector<math::XYZTLorentzVector>* fakeMuons_;
  std::vector<int>* fakeMuonTypes_;
  std::vector<double>* fakeMuonDR_;
  std::vector<double>* fakeMuonDPt_;
};

FakeMuonAnalyzer::FakeMuonAnalyzer(const edm::ParameterSet& pset)
{
  masses_ = new std::vector<double>();
  fakeMuons_ = new std::vector<math::XYZTLorentzVector>();
  fakeMuonTypes_ = new std::vector<int>();
  fakeMuonDR_ = new std::vector<double>();
  fakeMuonDPt_ = new std::vector<double>();

  muonLabel_ = pset.getParameter<edm::InputTag>("muon");
  kshortLabel_ = pset.getParameter<edm::InputTag>("kshort");

  std::vector<std::string> idCuts = pset.getParameter<std::vector<std::string> >("idCuts");
  for ( int i=0, n=idCuts.size(); i<n; ++i )
  {
    std::string idCut = idCuts[i];

    StringCutObjectSelector<reco::Muon, true>* idCutSel = new StringCutObjectSelector<reco::Muon, true>(idCut);
    idCuts_.push_back(idCutSel);
  }

  edm::Service<TFileService> fs;
  hEvent_ = fs->make<TH1F>("hEvent", "Event count", 5, 0, 5);
  hEvent_->GetXaxis()->SetBinLabel(1, "Total");
  hEvent_->GetXaxis()->SetBinLabel(2, "Muon");
  hEvent_->GetXaxis()->SetBinLabel(3, "Kshort");
  hEvent_->GetXaxis()->SetBinLabel(4, "Overlap");

  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");

  tree_->Branch("nMuon", &nMuon_, "nMuon/I");
  tree_->Branch("nKshort", &nKshort_, "nKshort/I");

  tree_->Branch("mass", "std::vector<double>", &masses_);
  tree_->Branch("fakeMuon", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &fakeMuons_);
  tree_->Branch("fakeMuonType", "std::vector<int>", &fakeMuonTypes_);
  tree_->Branch("fakeMuonDR", "std::vector<double>", &fakeMuonDR_);
  tree_->Branch("fakeMuonDPt", "std::vector<double>", &fakeMuonDPt_);

}

void FakeMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  masses_->clear();
  fakeMuons_->clear();
  fakeMuonTypes_->clear();
  fakeMuonDR_->clear();
  fakeMuonDPt_->clear();

  run_ = event.id().run();
  event_ = event.id().event();
  lumi_ = event.id().luminosityBlock();

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  edm::Handle<edm::View<reco::VertexCompositeCandidate> > kshortHandle;
  event.getByLabel(kshortLabel_, kshortHandle);

  nMuon_ = muonHandle->size();
  nKshort_ = kshortHandle->size();

  hEvent_->Fill(0);
  if ( nMuon_ != 0 ) hEvent_->Fill(1);
  if ( nKshort_ != 0 ) hEvent_->Fill(2);

  for ( int i=0; i<nMuon_; ++i )
  {
    const reco::Muon& iMuon = muonHandle->at(i);
    const math::XYZTLorentzVector muonLVec = iMuon.p4();
    int muonType = 0;
    for ( int j=0, nIdCuts=idCuts_.size(); j<nIdCuts; ++j )
    {
      if ( !(*idCuts_[j])(iMuon) ) continue;
      muonType |= 1<<j;
    }

    for ( int j=0; j<nKshort_; ++j )
    {
      const reco::CompositeCandidate& iKshort = kshortHandle->at(j);

      const reco::Candidate* p1 = iKshort.daughter(0);
      const reco::Candidate* p2 = iKshort.daughter(1);

      const reco::Candidate* pOverlap = 0;

      if ( iMuon.track() == p1->get<reco::TrackRef>() ) pOverlap = p1;
      else if ( iMuon.track() == p2->get<reco::TrackRef>() ) pOverlap = p2;

      if ( pOverlap )
      {
        //cout << "Muon: " << iMuon.track().id() << "\t" << iMuon.p4() << endl;
        //cout << "P1  : " << p1->get<reco::TrackRef>().id() << "\t" << p1->p4() << endl;
        //cout << "P2  : " << p2->get<reco::TrackRef>().id() << "\t" << p2->p4() << endl;
        masses_->push_back(iKshort.mass());
        fakeMuons_->push_back(muonLVec);
        fakeMuonTypes_->push_back(muonType);
        fakeMuonDR_->push_back(deltaR(*pOverlap, iMuon));
        fakeMuonDPt_->push_back(iMuon.pt() - pOverlap->pt());
      }
    }
  }

  if ( fakeMuons_->size() != 0 ) hEvent_->Fill(3);

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FakeMuonAnalyzer);
