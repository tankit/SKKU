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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
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

namespace muon
{
  reco::Muon::ArbitrationType arbitrationTypeFromString(const std::string& name)
  {
    if ( name == "NoArbitration" ) return reco::Muon::NoArbitration;
    else if ( name == "SegmentArbitration" ) return reco::Muon::SegmentArbitration;
    else if ( name == "SegmentAndTrackArbitration" ) return reco::Muon::SegmentAndTrackArbitration;
    else if ( name == "SegmentAndTrackArbitrationCleaned" ) return reco::Muon::SegmentAndTrackArbitrationCleaned;
    else if ( name == "RPCHitAndTrackArbitration" ) return reco::Muon::RPCHitAndTrackArbitration;
    else return reco::Muon::NoArbitration;
  }
}

class Histograms
{
public:
  Histograms(TFileDirectory* dir, const double massMin, const double massMax);
  void fill(const math::XYZTLorentzVector& lv, const double mass, const bool pass);
  
private:
  TH2F* hFrame_;
  std::map<std::string, TH1F*> hBinToMassMapPass_;
  std::map<std::string, TH1F*> hBinToMassMapFail_;

};

Histograms::Histograms(TFileDirectory* dir, const double massMin, const double massMax)
{
  const int nBinMass = int((massMax-massMin)*1000);

  const double ptBins[] = {0, 3, 5, 7, 10, 20, 30, 50, 150};
  const double etaBins[] = {-2.5, -2.1, -1.6, -1.2, -0.8, 0, 0.8, 1.2, 1.6, 2.1, 2.5};

  const int nPtBin = sizeof(ptBins)/sizeof(const double)-1;
  const int nEtaBin = sizeof(etaBins)/sizeof(const double)-1;

  hFrame_ = dir->make<TH2F>("hFrame", ";Transverse momentum p_{T} (GeV/c);Pseudorapidity #eta", nPtBin, ptBins, nEtaBin, etaBins);
  for ( int ptBin = 0; ptBin<nPtBin; ++ptBin )
  {
    for ( int etaBin = 0; etaBin<nEtaBin; ++etaBin )
    {
      TString binStr = Form("pt_bin%d__eta_bin%d", ptBin, etaBin);
      TFileDirectory binDir = dir->mkdir(binStr.Data());
      TString hTitle = Form(";Mass (GeV/c^{2});Entries per 1MeV/c^{2}");
      hBinToMassMapPass_[binStr.Data()] = binDir.make<TH1F>("hMass_Pass", "Pass"+hTitle, nBinMass, massMin, massMax);
      hBinToMassMapFail_[binStr.Data()] = binDir.make<TH1F>("hMass_Fail", "Fail"+hTitle, nBinMass, massMin, massMax);
    }
  }
}

void Histograms::fill(const math::XYZTLorentzVector& lv, const double mass, const bool pass)
{
  const double pt = lv.pt();
  const double eta = lv.eta();

  const int ptBin  = hFrame_->GetXaxis()->FindBin(pt) - 1;
  const int etaBin = hFrame_->GetYaxis()->FindBin(eta) - 1;

  TString binStr = Form("pt_bin%d__eta_bin%d", ptBin, etaBin);
  TH1F* h = 0;
  if ( pass ) h = hBinToMassMapPass_[binStr.Data()];
  else h = hBinToMassMapFail_[binStr.Data()];

  if ( h ) h->Fill(mass);
}

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
  std::vector<StringCutObjectSelector<reco::Muon, true>* > muonCuts_;
  std::vector<muon::SelectionType> muonSelectionTypes_;
  std::vector<reco::Muon::ArbitrationType> muonArbitrationTypes_;
  StringCutObjectSelector<reco::VertexCompositeCandidate, true>* vertexCut_;
  double maxDR_, maxDPt_;
  double massMin_, massMax_;

  bool doHist_;
  TH1F* hEvent_, * hNCandOverlap_;
  TTree* tree_;
  int run_, lumi_, event_;
  int nMuon_, nVertexCand_;
  int legId_;
  int muonCharge1_, muonCharge2_;
  int trackCharge1_, trackCharge2_;
  int pdgId_;
  double vertexMass_, vertexPt_, /*vertexL3D_, vertexL2D_,*/ vertexLxy_;
  //double deltaR_, deltaPt_;
  math::XYZTLorentzVector muon1_, muon2_;
  math::XYZTLorentzVector track1_, track2_;
  std::vector<int> muonIdResults1_, muonIdResults2_;
  //std::vector<int> muonStations_;
  std::vector<Histograms*> h1_, h2_;
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

  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  edm::Service<TFileService> fs;
  hEvent_ = fs->make<TH1F>("hEvent", "Event count", 5, 0, 5);
  hEvent_->GetXaxis()->SetBinLabel(1, "Total");
  hEvent_->GetXaxis()->SetBinLabel(2, "Muon");
  hEvent_->GetXaxis()->SetBinLabel(3, "Vertex");
  hEvent_->GetXaxis()->SetBinLabel(4, "Overlap");

  hNCandOverlap_ = fs->make<TH1F>("hNCandOverlap", "Number of vertex candidate with track sharing", 10, 0, 10);

  doHist_ = pset.getParameter<bool>("doHist");
  const bool doTree = pset.getParameter<bool>("doTree");
  tree_ = 0;
  if ( doTree )
  {
    tree_ = fs->make<TTree>("tree", "tree");
    tree_->Branch("run", &run_, "run/I");
    tree_->Branch("event", &event_, "event/I");
    tree_->Branch("lumi", &lumi_, "lumi/I");

    tree_->Branch("nMuon", &nMuon_, "nMuon/I");
    tree_->Branch("nVertex", &nVertexCand_, "nVertex/I");

    tree_->Branch("vertexMass", &vertexMass_, "vertexMass/D");
    tree_->Branch("pdgId", &pdgId_, "pdgId/I");
    tree_->Branch("vertexPt", &vertexPt_, "vertexPt/D");
    //tree_->Branch("vertexL3D", &vertexL3D_, "vertexL3D/D");
    //tree_->Branch("vertexL2D", &vertexL2D_, "vertexL2D/D");
    tree_->Branch("vertexLxy", &vertexLxy_, "vertexLxy/D");

    //tree_->Branch("deltaR" , &deltaR_ , "deltaR/D" );
    //tree_->Branch("deltaPt", &deltaPt_, "deltaPt/D");
    tree_->Branch("legId", &legId_, "legId/I");
    tree_->Branch("muonCharge1", &muonCharge1_, "muonCharge1/I");
    tree_->Branch("muonCharge2", &muonCharge2_, "muonCharge2/I");
    tree_->Branch("trackCharge1", &trackCharge1_, "trackCharge1/I");
    tree_->Branch("trackCharge2", &trackCharge2_, "trackCharge2/I");
    tree_->Branch("muon1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &muon1_ );
    tree_->Branch("muon2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &muon2_ );
    tree_->Branch("track1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &track1_);
    tree_->Branch("track2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &track2_);
  }

  edm::ParameterSet muonIds = pset.getParameter<edm::ParameterSet>("muonIds");
  std::vector<std::string> muonIdNames = muonIds.getParameterNames();
  muonIdResults1_.resize(muonIdNames.size());
  muonIdResults2_.resize(muonIdNames.size());
  //muonStations_.resize(muonIdNames.size());
  for ( int i=0, n=muonIdNames.size(); i<n; ++i )
  {
    const std::string& name = muonIdNames[i];
    edm::ParameterSet idCutSet = muonIds.getParameter<edm::ParameterSet>(name);
    const std::string cut = idCutSet.getParameter<std::string>("cut");
    const std::string idSelection = idCutSet.getParameter<std::string>("idSelection");
    const std::string arbitration = idCutSet.getParameter<std::string>("arbitration");

    muonCuts_.push_back(new StringCutObjectSelector<reco::Muon, true>(cut));
    muonSelectionTypes_.push_back(muon::selectionTypeFromString(idSelection));
    muonArbitrationTypes_.push_back(muon::arbitrationTypeFromString(arbitration));

    if ( doTree )
    {
      tree_->Branch(Form("muonId1_%s", name.c_str()), &(muonIdResults1_[0])+i, Form("muonId1_%s/I", name.c_str()));
      tree_->Branch(Form("muonId2_%s", name.c_str()), &(muonIdResults2_[0])+i, Form("muonId2_%s/I", name.c_str()));
      //tree_->Branch(Form("muStations_%s", name.c_str()), &(muonStations_[0])+i, Form("muStations_%s/I", name.c_str()));
    }

    if ( doHist_ )
    {
      TFileDirectory dir1 = fs->mkdir("hist/"+name+"/track1");
      h1_.push_back(new Histograms(&dir1, massMin_, massMax_));
      TFileDirectory dir2 = fs->mkdir("hist/"+name+"/track2");
      h2_.push_back(new Histograms(&dir2, massMin_, massMax_));
    }
  }
}

void FakeMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  run_ = event.id().run();
  event_ = event.id().event();
  lumi_ = event.id().luminosityBlock();

  for ( int i=0, n=muonIdResults1_.size(); i<n; ++i )
  {
    muonIdResults1_[i] = muonIdResults2_[i]  -1;
    //muonStations_[i] = -1;
  }

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
    vertexLxy_ = -1e9;
    //vertexL3D_ = vertexL2D_ = vertexLxy_ = -1e9;
    //deltaR_ = deltaPt_ = 1e-9;
    legId_ = 0;
    muonCharge1_ = muonCharge2_ = 0;
    trackCharge1_ = trackCharge2_ = 0;
    muon1_ = muon2_ = math::XYZTLorentzVector();
    track1_ = track2_ = math::XYZTLorentzVector();

    const reco::VertexCompositeCandidate& vertexCand = vertexCandHandle->at(iVertexCand);
    if ( !(*vertexCut_)(vertexCand) ) continue;

    pdgId_ = vertexCand.pdgId();
    vertexMass_ = vertexCand.mass();
    vertexPt_ = vertexCand.pt();

    if ( vertexMass_ < massMin_ or vertexMass_ > massMax_ ) continue;

    math::XYZPoint vertex = vertexCand.vertex();
    //vertexL3D_ = vertex.R();
    //vertexL2D_ = vertex.Rho();
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

    //const reco::Muon* matchedMuon = muon1 ? muon1 : muon2;
    //const reco::Candidate* matchedTrack = muon1 ? p1 : p2;
    // Set matched muon. 
    // legId = 0 for un-matched case.
    // There can be very rare case of double fake muons. 
    // Only muon1 will be stored in this case but we can separate this case with requiring legId != 3
    legId_ = 0;
    track1_ = p1->p4();
    trackCharge1_ = p1->charge();
    track2_ = p2->p4();
    trackCharge2_ = p2->charge();

    if ( muon1 ) 
    {
      legId_ |= 1;
      ++nFakeMuon; 
      muonCharge1_ = muon1->charge();
      muon1_ = muon1->p4();

      for ( int i=0, n=muonCuts_.size(); i<n; ++i )
      {
        muonIdResults1_[i] = true;
        muonIdResults1_[i] &= (*muonCuts_[i])(*muon1);
        muonIdResults1_[i] &= muon::isGoodMuon(*muon1, muonSelectionTypes_[i], muonArbitrationTypes_[i]);

      }
    }
    if ( muon2 ) 
    {
      legId_ |= 2;
      ++nFakeMuon; 
      muonCharge2_ = muon2->charge();
      muon2_ = muon2->p4();

      for ( int i=0, n=muonCuts_.size(); i<n; ++i )
      {
        muonIdResults2_[i] = true;
        muonIdResults2_[i] &= (*muonCuts_[i])(*muon2);
        muonIdResults2_[i] &= muon::isGoodMuon(*muon2, muonSelectionTypes_[i], muonArbitrationTypes_[i]);

      }
      
      //deltaR_ = deltaR(*matchedTrack, *matchedMuon);
      //deltaPt_ = matchedTrack->pt() - matchedMuon->pt();

      //muonStations_[i] = matchedMuon->numberOfMatchedStations(muonArbitrationTypes_[i]);
    }

    if ( doHist_ )
    {
      for ( int i=0, n=muonCuts_.size(); i<n; ++i )
      {
        h1_[i]->fill(track1_, vertexMass_, muonIdResults1_[i] == 1);
        h2_[i]->fill(track2_, vertexMass_, muonIdResults2_[i] == 1);
      }
    }
    if ( tree_ ) tree_->Fill();
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
