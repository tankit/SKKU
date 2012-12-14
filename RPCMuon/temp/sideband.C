#include<iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

void sideband()
{
//  TFile* f = TFile::Open("test_20121128.root");
  TFile* f = TFile::Open("result121212.root");
  TTree* tree = (TTree*)f->Get("fakeKshort/tree");
 

//  const char* cut_common = "abs(track.eta()) < 0.8";
  const char* cut_maxEta16 = "abs(track.eta()) < 1.6";
  const char* cut_barrel = "abs(track.eta()) <0.8";
  const char* cut_tran = "abs(track.eta()) >= 0.8 && abs(track.eta()) < 1.2";
  const char* cut_end = "abs(track.eta()) >= 1.2 && abs(track.eta()) < 1.6";

  double alpha = 1.0; 
  const char* cut_in = "abs(vertexMass-0.5) < 0.025";
  const char* cut_sb = "abs(vertexMass-0.5) > 0.040 && abs(vertexMass-0.5) < 0.065";

//  double alpha = 0.5;
//  const char* cut_in = "abs(vertexMass-0.5) < 0.025";
//  const char* cut_sb = "abs(vertexMass-0.5) > 0.050 && abs(vertexMass-0.5) < 0.1";

//  double alpha = 0.5;
//  const char* cut_in = "abs(vertexMass-0.5) < 0.025";
//  const char* cut_sb = "abs(vertexMass-0.5) > 0.040 && abs(vertexMass-0.5) < 0.09";


//  const char* cut_legId = "legId > -1";
  const char* cut_legId = "legId != 0";
 
  const char* cut_GLB = " muonId_globalMuon > 0";
  const char* cut_GLBMedium = " muonId_globalMuonMedium > 0";
  const char* cut_GLBTight = " muonId_globalMuonTight > 0";
  const char* cut_RPCMuLoose = " muonId_RPCMuLoose > 0";
  const char* cut_RPCMuMedium = "muonId_RPCMuMedium > 0";
  const char* cut_RPCMuPT = "muonId_RPCMuPromptTight > 0";
  const char* cut_TMLoose = "muonId_TMOneStationLoose > 0";
  const char* cut_TMTight = "muonId_TMOneStationTight > 0";
  const char* cut_TMTightPlus = "muonId_TMTwoStationTest > 0";
 
  TH1F* hMass_All = new TH1F("hMass_All", "Mass all;Mass", 100, 0.3, 0.7);
  TH1F* hMass_in = new TH1F("hMass_in" , "Inside;Mass", 100, 0.3, 0.7);
  TH1F* hMass_sb = new TH1F("hMass_sb", "sideBand;Mass", 100, 0.3, 0.7);

  TH1F* hPass_All = new TH1F("hPass_all", "Pass all;Pass", 100, 0.3, 0.7);
  TH1F* hPass_in = new TH1F("hPass_in", "inside;Pass", 100, 0.3, 0.7);
  TH1F* hPass_sb = new TH1F("hPass_sb", "sideband;Pass", 100, 0.3, 0.7); 

  TH1F* hMass_All_Barrel = new TH1F("hMass_All_Barrel", "Mass all;Mass", 100, 0.3, 0.7);
  TH1F* hMass_in_Barrel = new TH1F("hMass_in_Barrel" , "Inside;Mass", 100, 0.3, 0.7);
  TH1F* hMass_sb_Barrel = new TH1F("hMass_sb_Barrel", "sideBand;Mass", 100, 0.3, 0.7);

  TH1F* hPass_All_Barrel = new TH1F("hPass_all_Barrel", "Pass all;Pass", 100, 0.3, 0.7);
  TH1F* hPass_in_Barrel = new TH1F("hPass_in_Barrel", "inside;Pass", 100, 0.3, 0.7);
  TH1F* hPass_sb_Barrel = new TH1F("hPass_sb_Barrel", "sideband;Pass", 100, 0.3, 0.7);

  TH1F* hMass_All_Tran = new TH1F("hMass_All_Tran", "Mass all;Mass", 100, 0.3, 0.7);
  TH1F* hMass_in_Tran = new TH1F("hMass_in_Tran" , "Inside;Mass", 100, 0.3, 0.7);
  TH1F* hMass_sb_Tran = new TH1F("hMass_sb_Tran", "sideBand;Mass", 100, 0.3, 0.7);

  TH1F* hPass_All_Tran = new TH1F("hPass_all_Tran", "Pass all;Pass", 100, 0.3, 0.7);
  TH1F* hPass_in_Tran = new TH1F("hPass_in_Tran", "inside;Pass", 100, 0.3, 0.7);
  TH1F* hPass_sb_Tran = new TH1F("hPass_sb_Tran", "sideband;Pass", 100, 0.3, 0.7);

  TH1F* hMass_All_End = new TH1F("hMass_All_End", "Mass all;Mass", 100, 0.3, 0.7);
  TH1F* hMass_in_End = new TH1F("hMass_in_End" , "Inside;Mass", 100, 0.3, 0.7);
  TH1F* hMass_sb_End = new TH1F("hMass_sb_End", "sideBand;Mass", 100, 0.3, 0.7);

  TH1F* hPass_All_End = new TH1F("hPass_all_End", "Pass all;Pass", 100, 0.3, 0.7);
  TH1F* hPass_in_End = new TH1F("hPass_in_End", "inside;Pass", 100, 0.3, 0.7);
  TH1F* hPass_sb_End = new TH1F("hPass_sb_End", "sideband;Pass", 100, 0.3, 0.7);


/*
  //pt//
  double ptBins[15] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 15., 20., 25., 30.};
  TH1F* hPt_all_in = new TH1F("hPt_all_in", "In;P_{T}", 14, ptBins);
  TH1F* hPt_all_sb = new TH1F("hPt_all_sb", "In;P_{T}", 14, ptBins);
  TH1F* hPt_pass_in = new TH1F("hPt_pass_in", "In;P_{T}", 14, ptBins);
  TH1F* hPt_pass_sb = new TH1F("hPt_pass_sb", "In;P_{T}", 14, ptBins);
 
  TH1F* hPt_all = new TH1F("hPt_all", "All, Side band subtracted p_{T};P_{T}", 14, ptBins);
  TH1F* hPt_pass = new TH1F("hPt_pass", "Passing, Side band subtracted p_{T};P_{T}", 14, ptBins);
//  TH1F* hPt_ratio = new TH1F("hPt_ratio", "Fake ratio(p_{T}) : RPCMuLoose;p_{T}", 14, ptBins);
//  TH1F* hPt_ratio = new TH1F("hPt_ratio", "Fake ratio(p_{T}) : RPCMuMedium;p_{T}", 14, ptBins);
  TH1F* hPt_ratio = new TH1F("hPt_ratio", "Fake ratio(p_{T}) : RPCMuPromptTight;p_{T}", 14, ptBins);
//  TH1F* hPt_ratio = new TH1F("hPt_ratio", "Fake ratio(p_{T}) : GlobalMuonTight;p_{T}", 14, ptBins);
//  TH1F* hPt_ratio = new TH1F("hPt_ratio", "Fake ratio(p_{T}) : TMOneStaionTight;p_{T}", 14, ptBins);


  hPt_all_in->Sumw2();
  hPt_all_sb->Sumw2();
  hPt_pass_in->Sumw2();
  hPt_pass_sb->Sumw2();
  hPt_all->Sumw2();
  hPt_pass->Sumw2();
  hPt_ratio->Sumw2();
*/

  //p//
  double pBins[15] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 15., 20., 25., 30.};
  TH1F* hP_all_in = new TH1F("hP_all_in", "In;P", 14, pBins);
  TH1F* hP_all_sb = new TH1F("hP_all_sb", "In;P", 14, pBins);
  TH1F* hP_pass_in = new TH1F("hP_pass_in", "In;P", 14, pBins);
  TH1F* hP_pass_sb = new TH1F("hP_pass_sb", "In;P", 14, pBins);

  TH1F* hP_all = new TH1F("hP_all", "All, sideband subtracted p;p", 14, pBins);
  TH1F* hP_pass = new TH1F("hP_pass", "Passing, sideband subtracted p;p", 14, pBins);
//  TH1F* hP_ratio = new TH1F("hP_ratio", "Fake ratio(p) : RPCMuLoose;p", 14, pBins);
//  TH1F* hP_ratio = new TH1F("hP_ratio", "Fake ratio(p) : RPCMuMedium;p", 14, pBins);
  TH1F* hP_ratio = new TH1F("hP_ratio", "Fake ratio(p) : RPCMuPromptTight;p", 14, pBins);
//  TH1F* hP_ratio = new TH1F("hP_ratio", "Fake ratio(p) : GlobalMuonTight;p", 14, pBins);
//  TH1F* hP_ratio = new TH1F("hP_ratio", "Fake ratio(p) : TMOneStaionTightPlus;p", 14, pBins);

  hP_all_in->Sumw2();
  hP_all_sb->Sumw2();
  hP_pass_in->Sumw2();
  hP_pass_sb->Sumw2();
  hP_all->Sumw2();
  hP_pass->Sumw2();
  hP_ratio->Sumw2();

  //eta//
//  double etaBins[9] = {-2.5, -2.0, -1.5, -1.0, 0, 1.0, 1.5, 2.0, 2.5};
  double etaBins[7] = {-1.6, -1.2, -0.8, 0, 0.8, 1.2, 1.6};
  //double etaBins[3] = {0, 1.4, 2.5};
  TH1F* hEta_all_in = new TH1F("hEta_all_in", "In;#eta", 6, etaBins);
  TH1F* hEta_all_sb = new TH1F("hEta_all_sb", "In;#eta", 6, etaBins);
  TH1F* hEta_pass_in = new TH1F("hEta_pass_in", "In;#eta", 6, etaBins);
  TH1F* hEta_pass_sb = new TH1F("hEta_pass_sb", "In;#eta", 6, etaBins);

  TH1F* hEta_all = new TH1F("hEta_all", "All, Sideband subtracted #eta;#eta", 6, etaBins);
  TH1F* hEta_pass = new TH1F("hEta_pass", "Passing, Sideband subtracted #eta;#eta", 6, etaBins);
//  TH1F* hEta_ratio = new TH1F("hEta_ratio", "Fake ratio : RPCMuLoose;#eta", 6, etaBins);
//  TH1F* hEta_ratio = new TH1F("hEta_ratio", "Fake ratio : RPCMuMedium;#eta", 6, etaBins);
  TH1F* hEta_ratio = new TH1F("hEta_ratio", "Fake ratio : RPCMuPromptTight;#eta", 6, etaBins);
//  TH1F* hEta_ratio = new TH1F("hEta_ratio", "Fake ratio : GlobalMuonTight;#eta", 6, etaBins);
//  TH1F* hEta_ratio = new TH1F("hEta_ratio", "Fake ratio : TMOneStaionTightPlus;#eta", 6, etaBins);

  hEta_all_in->Sumw2();
  hEta_all_sb->Sumw2();
  hEta_pass_in->Sumw2();
  hEta_pass_sb->Sumw2();
  hEta_all->Sumw2();
  hEta_pass->Sumw2();
  hEta_ratio->Sumw2(); 

  
  //eta all 
  tree->Draw("vertexMass>>hMass_All", Form("(%s)", cut_maxEta16), "goff");
  tree->Draw("vertexMass>>hMass_in" , Form("(%s) && (%s)", cut_maxEta16, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb", Form("(%s) && (%s)", cut_maxEta16, cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_Barrel", Form("(%s)", cut_barrel), "goff");
  tree->Draw("vertexMass>>hMass_in_Barrel" , Form("(%s) && (%s)", cut_barrel, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_Barrel", Form("(%s) && (%s)", cut_barrel, cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_Tran", Form("(%s)", cut_tran), "goff");
  tree->Draw("vertexMass>>hMass_in_Tran" , Form("(%s) && (%s)", cut_tran, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_Tran", Form("(%s) && (%s)", cut_tran, cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_End", Form("(%s)", cut_end), "goff");
  tree->Draw("vertexMass>>hMass_in_End" , Form("(%s) && (%s)", cut_end, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_End", Form("(%s) && (%s)", cut_end, cut_sb), "goff");


  tree->Draw("vertexMass>>hPass_All", Form("(%s) && (%s) && (%s)", cut_maxEta16, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in", Form("(%s) && (%s) && (%s) && (%s)", cut_maxEta16, cut_RPCMuPT, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb", Form("(%s) && (%s) && (%s) && (%s)", cut_maxEta16, cut_RPCMuPT, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_Barrel", Form("(%s) && (%s) && (%s)", cut_barrel, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_Barrel", Form("(%s) && (%s) && (%s) && (%s)", cut_barrel, cut_RPCMuPT, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_Barrel", Form("(%s) && (%s) && (%s) && (%s)", cut_barrel, cut_RPCMuPT, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_Tran", Form("(%s) && (%s) && (%s)", cut_tran, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_Tran", Form("(%s) && (%s) && (%s) && (%s)", cut_tran, cut_RPCMuPT, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_Tran", Form("(%s) && (%s) && (%s) && (%s)", cut_tran, cut_RPCMuPT, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_End", Form("(%s) && (%s) && (%s)", cut_end, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_End", Form("(%s) && (%s) && (%s) && (%s)", cut_end, cut_RPCMuPT, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_End", Form("(%s) && (%s) && (%s) && (%s)", cut_end, cut_RPCMuPT, cut_legId, cut_sb), "goff");

/*
  tree->Draw("track.pt()>>hPt_all_in", Form("(%s)", cut_in), "goff");
  tree->Draw("track.pt()>>hPt_all_sb", Form("(%s)", cut_sb), "goff");
  tree->Draw("track.pt()>>hPt_pass_in", Form("(%s) && (%s) && (%s)", cut_in, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("track.pt()>>hPt_pass_sb", Form("(%s) && (%s) && (%s)", cut_sb, cut_RPCMuPT, cut_legId), "goff");
*/
  tree->Draw("track.P()>>hP_all_in", Form("(%s)", cut_in), "goff");
  tree->Draw("track.P()>>hP_all_sb", Form("(%s)", cut_sb), "goff");
  tree->Draw("track.P()>>hP_pass_in", Form("(%s) && (%s) && (%s)", cut_in, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("track.P()>>hP_pass_sb", Form("(%s) && (%s) && (%s)", cut_sb, cut_RPCMuPT, cut_legId), "goff");

  tree->Draw("track.eta()>>hEta_all_in", Form("(%s)", cut_in), "goff");
  tree->Draw("track.eta()>>hEta_all_sb", Form("(%s)", cut_sb), "goff");
  tree->Draw("track.eta()>>hEta_pass_in", Form("(%s) && (%s) && (%s)", cut_in, cut_RPCMuPT, cut_legId), "goff");
  tree->Draw("track.eta()>>hEta_pass_sb", Form("(%s) && (%s) && (%s)", cut_sb, cut_RPCMuPT, cut_legId), "goff");


/*  
  hPt_all->Add(hPt_all_in, hPt_all_sb, 1, -1);
  hPt_pass->Add(hPt_pass_in, hPt_pass_sb, 1, -1);
  hPt_ratio->Divide(hPt_pass, hPt_all, 1, 1, "B");
*/
  hP_all->Add(hP_all_in, hP_all_sb, 1, -1);
  hP_pass->Add(hP_pass_in, hP_pass_sb, 1, -1);
  hP_ratio->Divide(hP_pass, hP_all, 1, 1, "B");
//  hP_ratio->Divide(hP_pass, hP_all, 1, 1, "");
  
  hEta_all->Add(hEta_all_in, hEta_all_sb, 1, -1);
  hEta_pass->Add(hEta_pass_in, hEta_pass_sb, 1, -1);
  hEta_ratio->Divide(hEta_pass, hEta_all, 1, 1, "B");
//  hEta_ratio->Divide(hEta_pass, hEta_all, 1, 1, "");
 
  hMass_in->SetFillColor(kGreen);
  hMass_sb->SetFillColor(kBlue);

  //fake rate

  double sum1 = hMass_in->Integral();
  double sum2 = hMass_sb->Integral();
  double sum3 = hPass_in->Integral();
  double sum4 = hPass_sb->Integral();
  double sub1 = sum1 - (alpha*sum2);
  double sub2 = sum3 - (alpha*sum4);
  double R = sub2/sub1;

  double err1 = sqrt(((1-(2*R))*(sum3+sum4)+(R*R)*(sum1+sum2))/((sum1-sum2)*(sum1-sum2)));


  cout << "hMass_in:" << sum1 << ", hMass_sb:" << sum2 << ", hPass_in:" << sum3 << ", hPass_sb:" << sum4 <<endl;
  cout << "Fake Rate : " << R << endl;
  cout << "Error :" << err1 << endl;

  double sum_ba1 = hMass_in_Barrel->Integral();
  double sum_ba2 = hMass_sb_Barrel->Integral();
  double sum_ba3 = hPass_in_Barrel->Integral();
  double sum_ba4 = hPass_sb_Barrel->Integral();
  double sub_ba1 = sum_ba1 - (alpha*sum_ba2);
  double sub_ba2 = sum_ba3 - (alpha*sum_ba4);
  double R_ba = sub_ba2/sub_ba1;

  double err2 = sqrt(((1-(2*R_ba))*(sum_ba3+sum_ba4)+(R_ba*R_ba)*(sum_ba1+sum_ba2))/((sum_ba1-sum_ba2)*(sum_ba1-sum_ba2)));


  cout << "hMass_in_barrel :" << sum_ba1 << ", hMass_sb_barrel:" << sum_ba2 << ", hPass_in_barrel:" << sum_ba3 << ", hPass_sb_barrel:" << sum_ba4 << endl;
  cout << "Fake Rate(Barrel): " << R_ba << endl;
  cout << "Error : " << err2 << endl;

  double sum_tr1 = hMass_in_Tran->Integral();
  double sum_tr2 = hMass_sb_Tran->Integral();
  double sum_tr3 = hPass_in_Tran->Integral();
  double sum_tr4 = hPass_sb_Tran->Integral();
  double sub_tr1 = sum_tr1 - (alpha*sum_tr2);
  double sub_tr2 = sum_tr3 - (alpha*sum_tr4);
  double R_tr = sub_tr2/sub_tr1;

  double err3 = sqrt(((1-(2*R_tr))*(sum_tr3+sum_tr4)+(R_tr*R_tr)*(sum_tr1+sum_tr2))/((sum_tr1-sum_tr2)*(sum_tr1-sum_tr2)));


  cout << "hMass_in_tran :" << sum_tr1 << ", hMass_sb_tran:" << sum_tr2 << ", hPass_in_tran:" << sum_tr3 << ", hPass_sb_tran:" << sum_tr4 << endl;
  cout << "Fake Rate(Transition): " << R_tr << endl;
  cout << "Error : " << err3 << endl;

  double sum_end1 = hMass_in_End->Integral();
  double sum_end2 = hMass_sb_End->Integral();
  double sum_end3 = hPass_in_End->Integral();
  double sum_end4 = hPass_sb_End->Integral();
  double sub_end1 = sum_end1 - (alpha*sum_end2);
  double sub_end2 = sum_end3 - (alpha*sum_end4);
  double R_end = sub_end2/sub_end1;

  double err4 = sqrt(((1-(2*R_end))*(sum_end3+sum_end4)+(R_end*R_end)*(sum_end1+sum_end2))/((sum_end1-sum_end2)*(sum_end1-sum_end2)));


  cout << "hMass_in_end :" << sum_end1 << ", hMass_sb_end:" << sum_end2 << ", hPass_in_end:" << sum_end3 << ", hPass_sb_end:" << sum_end4 << endl;
  cout << "Fake Rate(Endcap): " << R_end << endl;
  cout << "Error : " << err4 << endl;

 //Draw 
  TCanvas* cMass = new TCanvas("cMass", "cMass", 500, 500);
  hMass_All->Draw();
  hMass_in->Draw("same");
  hMass_sb->Draw("same");
 
//  TCanvas* cPt = new TCanvas("cPt", "cPt", 500, 500);
//  hPt_ratio->Draw();

  TCanvas* cP = new TCanvas("cP", "cP", 500, 500);
  hP_ratio->Draw();

  TCanvas* cEta = new TCanvas("cEta", "cEta", 500, 500);
  hEta_ratio->Draw("e");
  //cEta->Print("_eta.png");

}
