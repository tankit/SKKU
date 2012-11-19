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

//const char* sample = "Summer12_TTJets";
const char* sample = "Run2012B_MuHad";
//const char* sample = "Fall11_TTJets";
//const char* sample = "Run2011A_SingleMu";

const int cut_minNJet = 4;
const double cut_minJetPt = 35;
const double cut_maxJetEta = 2.5;
const double cut_minLeadJetPt = 45;

const int cut_minNBjet = 2;
const double cut_bTag = 0.679;

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
  TFile* file1 = TFile::Open(Form("ntuple/result_%s.root", sample));
  TFile* file2 = TFile::Open(Form("ntuple/result_%s.root", sample));
  TTree* tree1 = (TTree*)file1->Get("event/tree");
  TTree* tree2 = (TTree*)file2->Get("event/tree");

  EventTopology event1, event2;

  setBranch(tree1, event1);
  setBranch(tree2, event2);

  TFile* outFile = TFile::Open(Form("hist_%s_pt%.0f_nj%d_nb%d.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet), "RECREATE");
  TDirectory* dirSEvt = outFile->mkdir("SEvt");
  TDirectory* dirBEvt = outFile->mkdir("BEvt");

  dirSEvt->cd();
  const double binShift = 2.5; //0;
  const double mJJMax = 500+binShift;
  const double mJJBMax = 1000+binShift;
  const int nbin = 100;
  const TString axisTitleMw = ";Dijet mass (GeV/c^{2});Events per 5GeV/c^{2}";
  const TString axisTitleMt = ";Three jet mass (GeV/c^{2});Events per 10GeV/c^{2}";

  TH1F* hSEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* hSEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, binShift, mJJBMax);
  TH1F* hSEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hSEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hSEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hSEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, binShift, mJJBMax);
  TH1F* hSEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, binShift, mJJBMax);

  dirBEvt->cd();
  TH1F* hBEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* hBEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, binShift, mJJBMax);
  TH1F* hBEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hBEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hBEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, binShift, mJJBMax);
  //TH1F* hBEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, binShift, mJJBMax);
  TH1F* hBEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, binShift, mJJBMax);
  
  int nBiEvent = 0, nPassedBiEvent = 0; // Variables to restore overlap removal scale
  for ( int iEvent1=0, nEvent=tree1->GetEntries(); iEvent1<nEvent; ++iEvent1 )
  {
    tree1->GetEntry(iEvent1);
    const LorentzVector& lepton1 = *event1.lepton;

    if ( event1.jets->size() < cut_minNJet ) continue;
    if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
    int nBjet1 = 0;
    for ( int j1=0, nj=event1.bTags->size(); j1<nj; ++j1 )
    {
      if ( event1.bTags->at(j1) > cut_bTag ) ++nBjet1;
    }
    if ( nBjet1 < cut_minNBjet ) continue;

    // Make SEvt event jet combinations
    for ( int j1=0, nj=event1.jets->size(); j1<nj; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);
      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;

      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      //if ( mcBit1&3 ) continue;
      if ( bTag1 > cut_bTag ) continue;
      for ( int j2=j1+1; j2<nj; ++j2 )
      {
        const LorentzVector jet2 = event1.jets->at(j2);
        const double bTag2 = event1.bTags->at(j2);
        const int mcBit2 = event1.jetMCBits->at(j2);

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;

        if ( DeltaR(jet1, jet2) < 0.4 ) continue;

        LorentzVector jj = event1.jets->at(j1)+event1.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hSEvt_Mw->Fill(mJJ);

        if ( mcBit1&4 and mcBit2&4 ) hSEvt_Mw_JJ->Fill(mJJ);
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hSEvt_Mw_HBJ->Fill(mJJ);
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) hSEvt_Mw_LBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hSEvt_Mw_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hSEvt_Mw_KK->Fill(mJJ);
        else hSEvt_Mw_JK->Fill(mJJ);

        for ( int j3=0; j3<nj; ++j3 )
        {
          const LorentzVector jet3 = event1.jets->at(j3);
          const double bTag3 = event1.bTags->at(j3);
          const int mcBit3 = event1.jetMCBits->at(j3);

          //if ( DeltaR(jet3, lepton1) < 0.3 ) continue;
          if ( jet3.pt() < cut_minJetPt or abs(jet3.eta()) > cut_maxJetEta ) continue;
          if ( bTag3 <= cut_bTag ) continue;
          if ( DeltaR(jet1, jet3) < 0.4 ) continue;
          if ( DeltaR(jet2, jet3) < 0.4 ) continue;

          LorentzVector jjb = jj + event1.jets->at(j3);
          const double mJJB = jjb.mass();

          hSEvt_Mt->Fill(mJJB);

          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hSEvt_Mt_JJHB->Fill(mJJB);
          //else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
          //          (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
          //          (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hSEvt_Mt_JJLB->Fill(mJJB);
          //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hSEvt_Mt_KKX->Fill(mJJB);
          //else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hSEvt_Mt_KXY->Fill(mJJB);
          else hSEvt_Mt_XYZ->Fill(mJJB);
        }
      }
      //break;

    }

    // Make Bi event jet combinations
    for ( int iEvent2 = iEvent1+1; ;++iEvent2 )
    {
      if ( iEvent2 == nEvent ) iEvent2 = 0;
      tree2->GetEntry(iEvent2);

      const LorentzVector& lepton2 = *event2.lepton;

      if ( event2.jets->size() < cut_minNJet ) continue;
      if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
      int nBjet2 = 0;
      for ( int j2=0, nj=event2.bTags->size(); j2<nj; ++j2 )
      {
        if ( event2.bTags->at(j2) > cut_bTag ) ++nBjet2;
      }
      if ( nBjet2 < cut_minNBjet ) continue;

      break;
    }

    for ( int j1=0, nj1=event1.jets->size(); j1<nj1; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);

      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;
      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      if ( bTag1 > cut_bTag ) continue;
      //if ( mcBit1&3 ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;

      for ( int j2=0, nj2=event2.jets->size(); j2<nj2; ++j2 )
      {
        const LorentzVector jet2 = event2.jets->at(j2);
        const double bTag2 = event2.bTags->at(j2);
        const int mcBit2 = event2.jetMCBits->at(j2);

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;

        ++nBiEvent;
        if ( DeltaR(jet1, jet2) < 0.4 ) continue;
        ++nPassedBiEvent;

        LorentzVector jj = event1.jets->at(j1)+event2.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hBEvt_Mw->Fill(mJJ);

        if ( mcBit1&4 and mcBit2&4 ) hBEvt_Mw_JJ->Fill(mJJ);
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hBEvt_Mw_HBJ->Fill(mJJ);
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) hBEvt_Mw_LBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hBEvt_Mw_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hBEvt_Mw_KK->Fill(mJJ);
        else hBEvt_Mw_JK->Fill(mJJ);

        for ( int j3=0; j3<nj2; ++j3 )
        {
          const LorentzVector jet3 = event2.jets->at(j3);
          const double bTag3 = event2.bTags->at(j3);
          const int mcBit3 = event2.jetMCBits->at(j3);

          //if ( DeltaR(jet3, lepton1) < 0.3 ) continue;
          if ( jet3.pt() < cut_minJetPt or abs(jet3.eta()) > cut_maxJetEta ) continue;
          if ( bTag3 <= cut_bTag ) continue;
          if ( DeltaR(jet1, jet3) < 0.4 ) continue;
          if ( DeltaR(jet2, jet3) < 0.4 ) continue;

          LorentzVector jjb = jj + event2.jets->at(j3);
          const double mJJB = jjb.mass();

          hBEvt_Mt->Fill(mJJB);

          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB);
          //else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
          //          (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
          //          (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hBEvt_Mt_JJLB->Fill(mJJB);
          //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hBEvt_Mt_KKX->Fill(mJJB);
          //`else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hBEvt_Mt_KXY->Fill(mJJB);
          else hBEvt_Mt_XYZ->Fill(mJJB);
         
        }
      }
      //break;

    }


    //cout << event1.lepton->pt() << ' ' << event2.lepton->pt() << endl;
    //cout << event1.charge << ' ' << event2.charge << endl;
    //cout << event1.bTags->size() << ' ' << event2.bTags->size() << endl;
  }

  // Restore scale factor by skipping overlapping jets
  //hBEvt_Mw->Scale(0.5*nBiEvent/nPassedBiEvent);
  //hBEvt_Mw->Scale(0.5);
  //hBEvt_Mt->Scale(1./3);

  outFile->Write();

}
