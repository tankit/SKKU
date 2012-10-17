#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/VectorUtil_Cint.h"

using namespace std;
using namespace ROOT::Math::VectorUtil;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

const int cut_minNJet = 3;
const double cut_minJetPt = 40;
const double cut_minLeadJetPt = 45;

const int cut_minNBjet = 1;
const double cut_bTag = 0.679;
//const double cut_bTag = 0.2; // Arbitrary tight b-veto cut

struct EventTopology
{
  int run, lumi, event;
  LorentzVector* lepton;
  LorentzVector* met;
  std::vector<LorentzVector>* jets;
  int charge;
  std::vector<double>* bTags;
  std::vector<int>* jetMCBits;

  EventTopology()
  {
    lepton = new LorentzVector;
    met = new LorentzVector;
    jets = new std::vector<LorentzVector>();
    bTags = new std::vector<double>();
    jetMCBits = new std::vector<int>();
  };
};

void setBranch(TTree* tree, EventTopology& e)
{
  tree->SetBranchAddress("run"  , &e.run  );
  tree->SetBranchAddress("lumi" , &e.lumi );
  tree->SetBranchAddress("event", &e.event);

  tree->SetBranchAddress("lepton", &e.lepton);
  tree->SetBranchAddress("met", &e.met);
  tree->SetBranchAddress("jets", &e.jets);

  tree->SetBranchAddress("charge", &e.charge);
  tree->SetBranchAddress("bTags", &e.bTags);
  tree->SetBranchAddress("jetMCBits", &e.jetMCBits);

}

void mix()
{
  TFile* file1 = TFile::Open("/users/jwseo/work/CMS/CMSSW_5_2_6/src/SKKU/TTbarLeptonJet/test/result_1004.root");
  TFile* file2 = TFile::Open("/users/jwseo/work/CMS/CMSSW_5_2_6/src/SKKU/TTbarLeptonJet/test/result_1004.root");
  TTree* tree1 = (TTree*)file1->Get("event/tree");
  TTree* tree2 = (TTree*)file2->Get("event/tree");

  EventTopology event1, event2;

  setBranch(tree1, event1);
  setBranch(tree2, event2);

  TFile* outFile = TFile::Open("hist.root", "RECREATE");
  TDirectory* dirSevt = outFile->mkdir("SEvt");
  TDirectory* dirBevt = outFile->mkdir("BEvt");

  dirSevt->cd();
  const double mJJMax = 500+2.5;
  const int nbin = 100;
  TH1F* hSevt_M_JJ = new TH1F("hM_JJ", "Dijet mass;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_JJ = new TH1F("hMMC_JJ", "J+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_JK = new TH1F("hMMC_JK", "J+K;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_KK = new TH1F("hMMC_KK", "K+K;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_HBJ = new TH1F("hMMC_HBJ", "HadB+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_LBJ = new TH1F("hMMC_LBJ", "LepB+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hSevt_MMC_BX = new TH1F("hMMC_BX", "B+X;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  //TH1F* hSevt_MMC_JB = new TH1F("hMMC_JB", "J+b;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  //TH1F* hSevt_
  TH1F* hSevt_M_JJB = new TH1F("hM_JJB", "Three jet mass;M(jjB)", nbin, 2.5, 1000);

  dirBevt->cd();
  TH1F* hBevt_M_JJ = new TH1F("hM_JJ", "Dijet mass;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_JJ = new TH1F("hMMC_JJ", "J+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_JK = new TH1F("hMMC_JK", "J+K;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_KK = new TH1F("hMMC_KK", "K+K;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_HBJ = new TH1F("hMMC_HBJ", "HadB+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_LBJ = new TH1F("hMMC_LBJ", "LepB+J;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);
  TH1F* hBevt_MMC_BX = new TH1F("hMMC_BX", "B+X;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", nbin, 2.5, mJJMax);

  TH1F* hBevt_M_JJB = new TH1F("hM_JJB", "Three jet mass;M(jjB)", nbin, 2.5, 1000);
  
  int nBiEvent = 0, nPassedBiEvent = 0; // Variables to restore overlap removal scale
  for ( int iEvent=0, nEvent=tree1->GetEntries(); iEvent<nEvent; ++iEvent )
  {
    tree1->GetEntry(iEvent);

    if ( event1.jets->size() < cut_minNJet ) continue;
    if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
    int nBjet1 = 0;
    for ( int j1=0, nj=event1.bTags->size(); j1<nj; ++j1 )
    {
      if ( event1.bTags->at(j1) > cut_bTag ) ++nBjet1;
    }
    if ( nBjet1 < cut_minNBjet ) continue;

    // Make Sevt event jet combinations
    for ( int j1=0, nj=event1.jets->size(); j1<nj; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);

      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > 2.5 ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      //if ( mcBit1&3 ) continue;
      if ( bTag1 > cut_bTag ) continue;
      for ( int j2=j1+1; j2<nj; ++j2 )
      {
        const LorentzVector jet2 = event1.jets->at(j2);
        const double bTag2 = event1.bTags->at(j2);
        const int mcBit2 = event1.jetMCBits->at(j2);

        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > 2.5 ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;

        if ( DeltaR(jet1, jet2) < 0.4 ) continue;

        LorentzVector jj = event1.jets->at(j1)+event1.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hSevt_M_JJ->Fill(mJJ);

        if ( (mcBit1&1 and mcBit2 != 0) or (mcBit2&1 and mcBit1 != 0) ) cout << "Leptonic b " << mcBit1 << " " << jet1.eta() << "," << jet1.phi() << " " << mcBit2 << " " << jet2.eta() << "," << jet2.phi() << "\n";

        if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) hSevt_MMC_LBJ->Fill(mJJ);
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hSevt_MMC_HBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hSevt_MMC_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hSevt_MMC_KK->Fill(mJJ);
        else if ( mcBit1&4 and mcBit2&4 ) hSevt_MMC_JJ->Fill(mJJ);
        else hSevt_MMC_JK->Fill(mJJ);

/*
        for ( int j3=j1+2; j3<nj; ++j3 )
        {
          const LorentzVector jet3 = event1.jets->at(j3);
          const double bTag3 = event1.bTags->at(j3);
          if ( bTag3 <= cut_bTag ) continue;
          //if ( (event1.jetMCBits->at(j3) & (1+2)) == 0 ) continue;
          //if ( (event1.jetMCBits->at(j3) & 4) != 0 ) continue;
          //if ( DeltaR(jet1, jet3) < 0.4 ) continue;
          //if ( DeltaR(jet2, jet3) < 0.4 ) continue;
          LorentzVector jjb = jj + event1.jets->at(j3);

          hSevt_M_JJB->Fill(jjb.mass());
        }
*/
      }
      //break;

    }

    // Make Bi event jet combinations
    tree2->GetEntry((iEvent+1)%nEvent);

    if ( event2.jets->size() < cut_minNJet ) continue;
    if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
    int nBjet2 = 0;
    for ( int j2=0, nj=event2.bTags->size(); j2<nj; ++j2 )
    {
      if ( event2.bTags->at(j2) > cut_bTag ) ++nBjet2;
    }
    if ( nBjet2 < cut_minNBjet ) continue;

    for ( int j1=0, nj1=event1.jets->size(); j1<nj1; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);

      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > 2.5 ) continue;
      if ( bTag1 > cut_bTag ) continue;
      //if ( mcBit1&3 ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;

      for ( int j2=0, nj2=event2.jets->size(); j2<nj2; ++j2 )
      {
        const LorentzVector jet2 = event2.jets->at(j2);
        const double bTag2 = event2.bTags->at(j2);
        const int mcBit2 = event2.jetMCBits->at(j2);

        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > 2.5 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;

        ++nBiEvent;
        if ( DeltaR(jet1, jet2) < 0.4 ) continue;
        ++nPassedBiEvent;

        LorentzVector jj = event1.jets->at(j1)+event2.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hBevt_M_JJ->Fill(mJJ);

        if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hBevt_MMC_HBJ->Fill(mJJ);
        else if ( (mcBit1==1 and mcBit2==4) or (mcBit1==4 and mcBit2==1) ) hBevt_MMC_LBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hBevt_MMC_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hBevt_MMC_KK->Fill(mJJ);
        else if ( mcBit1==4 and mcBit2==4 ) hBevt_MMC_JJ->Fill(mJJ);
        else hBevt_MMC_JK->Fill(mJJ);

        /*
        for ( int j3=j2+1; j3<nj2; ++j3 )
        {
          const LorentzVector jet3 = event2.jets->at(j3);
          const double bTag3 = event2.bTags->at(j3);
          //if ( bTag3 <= cut_bTag ) continue;
          //if ( (event2.jetMCBits->at(j3) & (1+2)) == 0 ) continue;
          //if ( (event2.jetMCBits->at(j3) & 4) != 0 ) continue;
          //if ( DeltaR(jet1, jet3) < 0.4 ) continue;
          //if ( DeltaR(jet2, jet3) < 0.4 ) continue;

          LorentzVector jjb = jj + event2.jets->at(j3);
          hBevt_M_JJB->Fill(jjb.mass());//, 1./4);
        }
        */
      }
      //break;

    }


    //cout << event1.lepton->pt() << ' ' << event2.lepton->pt() << endl;
    //cout << event1.charge << ' ' << event2.charge << endl;
    //cout << event1.bTags->size() << ' ' << event2.bTags->size() << endl;
  }

  // Restore scale factor by skipping overlapping jets
  //hBevt_M_JJ->Scale(0.5*nBiEvent/nPassedBiEvent);
  //hBevt_M_JJ->Scale(0.5);
  //hBevt_M_JJB->Scale(1./3);

  outFile->Write();

}
