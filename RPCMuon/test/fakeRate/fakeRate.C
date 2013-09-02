#include<iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
  
using namespace std; 

void fakeRate()
{ 
  TFile* f = TFile::Open("test121210.root");
//  TFile* f = TFile::Open("result_20121130.root");
  TTree* tree = (TTree*)f->Get("fakeKshort/tree");

//  const char* cut_common = "abs(track.eta()) < 0.8";
  const char* cut_barrel = "abs(track.eta()) <0.8";
  const char* cut_tran = "abs(track.eta()) >= 0.8 && abs(track.eta()) < 1.2";
  const char* cut_end = "abs(track.eta()) >= 1.2 && abs(track.eta()) < 1.6";
  const char* cut_in = "abs(vertexMass-0.5) < 0.025";
  const char* cut_sb = "abs(vertexMass-0.5) > 0.040 && abs(vertexMass-0.5) < 0.065";

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

  tree->Draw("vertexMass>>hMass_All", "", "goff");
  tree->Draw("vertexMass>>hMass_in" , Form("(%s)", cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb", Form("(%s)", cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_Barrel", Form("(%s)", cut_barrel), "goff");
  tree->Draw("vertexMass>>hMass_in_Barrel" , Form("(%s) && (%s)", cut_barrel, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_Barrel", Form("(%s) && (%s)", cut_barrel, cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_Tran", Form("(%s)", cut_tran), "goff");
  tree->Draw("vertexMass>>hMass_in_Tran" , Form("(%s) && (%s)", cut_tran, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_Tran", Form("(%s) && (%s)", cut_tran, cut_sb), "goff");

  tree->Draw("vertexMass>>hMass_All_End", Form("(%s)", cut_end), "goff");
  tree->Draw("vertexMass>>hMass_in_End" , Form("(%s) && (%s)", cut_end, cut_in), "goff");
  tree->Draw("vertexMass>>hMass_sb_End", Form("(%s) && (%s)", cut_end, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All", Form("(%s) && (%s)",cut_TMLoose, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in", Form("(%s) && (%s) && (%s)", cut_TMLoose, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb", Form("(%s) && (%s) && (%s)", cut_TMLoose, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_Barrel", Form("(%s) && (%s) && (%s)", cut_barrel, cut_TMLoose, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_Barrel", Form("(%s) && (%s) && (%s) && (%s)", cut_barrel, cut_TMLoose, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_Barrel", Form("(%s) && (%s) && (%s) && (%s)", cut_barrel, cut_TMLoose, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_Tran", Form("(%s) && (%s) && (%s)", cut_tran, cut_TMLoose, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_Tran", Form("(%s) && (%s) && (%s) && (%s)", cut_tran, cut_TMLoose, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_Tran", Form("(%s) && (%s) && (%s) && (%s)", cut_tran, cut_TMLoose, cut_legId, cut_sb), "goff");

  tree->Draw("vertexMass>>hPass_All_End", Form("(%s) && (%s) && (%s)", cut_end, cut_TMLoose, cut_legId), "goff");
  tree->Draw("vertexMass>>hPass_in_End", Form("(%s) && (%s) && (%s) && (%s)", cut_end, cut_TMLoose, cut_legId, cut_in), "goff");
  tree->Draw("vertexMass>>hPass_sb_End", Form("(%s) && (%s) && (%s) && (%s)", cut_end, cut_TMLoose, cut_legId, cut_sb), "goff");

  //fake rate
  double sum1 = hMass_in->Integral();
  double sum2 = hMass_sb->Integral();
  double sum3 = hPass_in->Integral();
  double sum4 = hPass_sb->Integral();
  double sub1 = sum1 - sum2;
  double sub2 = sum3 - sum4;
  double R = sub2/sub1;

  double err1 = sqrt(((1-(2*R))*(sum3+sum4)+(R*R)*(sum1+sum2))/((sum1-sum2)*(sum1-sum2)));

  cout << "hMass_in:" << sum1 << "hMass_sb:" << sum2 << "hPass_in:" << sum3 << "hPass_sb:" << sum4 <<endl;
  cout << "Fake Rate : " << R << endl;
  cout << "Error :" << err1 << endl;

  double sum_ba1 = hMass_in_Barrel->Integral();
  double sum_ba2 = hMass_sb_Barrel->Integral();
  double sum_ba3 = hPass_in_Barrel->Integral();
  double sum_ba4 = hPass_sb_Barrel->Integral();
  double sub_ba1 = sum_ba1 - sum_ba2;
  double sub_ba2 = sum_ba3 - sum_ba4;
  double R_ba = sub_ba2/sub_ba1;

  double err2 = sqrt(((1-(2*R_ba))*(sum_ba3+sum_ba4)+(R_ba*R_ba)*(sum_ba1+sum_ba2))/((sum_ba1-sum_ba2)*(sum_ba1-sum_ba2)));

  cout << "hMass_in_barrel :" << sum_ba1 << "hMass_sb_barrel:" << sum_ba2 << "hPass_in_barrel:" << sum_ba3 << "hPass_sb_barrel:" << sum_ba4 << endl;
  cout << "Fake Rate(Barrel): " << R_ba << endl;
  cout << "Error : " << err2 << endl; 

  double sum_tr1 = hMass_in_Tran->Integral();
  double sum_tr2 = hMass_sb_Tran->Integral();
  double sum_tr3 = hPass_in_Tran->Integral();
  double sum_tr4 = hPass_sb_Tran->Integral();
  double sub_tr1 = sum_tr1 - sum_tr2;
  double sub_tr2 = sum_tr3 - sum_tr4;
  double R_tr = sub_tr2/sub_tr1;

  double err3 = sqrt(((1-(2*R_tr))*(sum_tr3+sum_tr4)+(R_tr*R_tr)*(sum_tr1+sum_tr2))/((sum_tr1-sum_tr2)*(sum_tr1-sum_tr2)));

  cout << "hMass_in_tran :" << sum_tr1 << "hMass_sb_tran:" << sum_tr2 << "hPass_in_tran:" << sum_tr3 << "hPass_sb_tran:" << sum_tr4 << endl;
  cout << "Fake Rate(Transition): " << R_tr << endl;
  cout << "Error : " << err3 << endl;

  double sum_end1 = hMass_in_End->Integral();
  double sum_end2 = hMass_sb_End->Integral();
  double sum_end3 = hPass_in_End->Integral();
  double sum_end4 = hPass_sb_End->Integral();
  double sub_end1 = sum_end1 - sum_end2;
  double sub_end2 = sum_end3 - sum_end4;
  double R_end = sub_end2/sub_end1;

  double err4 = sqrt(((1-(2*R_end))*(sum_end3+sum_end4)+(R_end*R_end)*(sum_end1+sum_end2))/((sum_end1-sum_end2)*(sum_end1-sum_end2)));

  cout << "hMass_in_end :" << sum_end1 << "hMass_sb_end:" << sum_end2 << "hPass_in_end:" << sum_end3 << "hPass_sb_end:" << sum_end4 << endl;
  cout << "Fake Rate(Endcap): " << R_end << endl;
  cout << "Error : " << err4 << endl;

}
