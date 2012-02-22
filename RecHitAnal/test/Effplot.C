#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <map>
#include <fstream.h>
////#include "/afs/cern.ch/user/m/mskim/public/styleTnP.h"

void Effplot(TString var="eta", Float_t hmin = 0.0){

  //gROOT->LoadMacro("/afs/cern.ch/user/m/mskim/public/tdrStyle.C");
  //setTDRStyle();

//  using namespace std;
//  ofstream fout; 
//  fout.open("Efficiency_Table.txt"); 
  cout << "Output of data" <<endl;

  TString ytitle = "Efficiency";

  //TString var = "muon_pt_eta";
  //TString var1 = "eta";
  //TString var2 = "pt_bin0";
  //TString varDir = var1+"_PLOT_"+var2

  TString varDir = var+"_PLOT";

  double Xmin = -1.8, Xmax = 1.8;
  Xmin = -2.0, Xmax = 2.0;

  TString FileIso="./data_eff.root", FileIso_mc="./mc_eff.root";

  TString DirIso="muonEffs/"+var+"/fit_eff_plots";

  TString xtitle = "Probe #eta";  TString htitle = "A RooPlot of Probe #eta";
  if(var=="pt") {
    Xmin = 20, Xmax = 100;
    xtitle = "Probe P_{t} (GeV/c)";
    htitle = "A RooPlot of Probe P_{t}";
  }

  cout << "///////////////////// Isolation (Data) ///////////////////////////" << endl;
  TFile * f_Iso = new TFile(Form("%s",FileIso.Data()));
  f_Iso->cd(Form("%s",DirIso.Data()));
    
  TGraphAsymmErrors * gr_Iso = new TGraphAsymmErrors();
  TCanvas* c1 = (TCanvas*) gDirectory->FindKey(varDir)->ReadObj();

  TString obj1="";
  obj1 = "hxy_fit_eff"; //hard-coded: eta_PLOT_pt_bin0->GetListOfPrimitives()->At(2)->GetName()
  RooHist* h1 = (RooHist*) c1->FindObject(obj1);
  int nbin = h1->GetMaxSize(); cout << nbin << endl;
  for(int j=0 ; j < nbin ; j++){
    double x;
    double y;
    double xerrhi = h1->GetErrorXhigh(j);
    double xerrlo = h1->GetErrorXlow(j);
    double yerrhi = h1->GetErrorYhigh(j);
    double yerrlo = h1->GetErrorYlow(j);
    int ibin =  h1->GetPoint(j,x,y);
    cout << "[" << x-xerrlo  << "," << x+xerrhi << "] "  << " eff (" << ibin << ") = " << y << " (+" << yerrhi << " -" << yerrlo << ")" << endl;
    gr_Iso->SetPoint(j,x,y);
    gr_Iso->SetPointError(j,xerrlo,xerrhi,yerrlo,yerrhi);          
  }
  
  //////////////////////////////////////////////////////////////////////////
  
  cout << "///////////////////// Isolation (MC) ///////////////////////////" << endl;
  TFile * f_Iso_mc = new TFile(Form("%s",FileIso_mc.Data()));
  f_Iso_mc->cd(Form("%s",DirIso.Data()));
    
  TGraphAsymmErrors * gr_Isomc = new TGraphAsymmErrors();
  TCanvas* c1 = (TCanvas*) gDirectory->FindKey(varDir)->ReadObj();
  TString obj1="";
  obj1 = "hxy_fit_eff";
  RooHist* h1 = (RooHist*) c1->FindObject(obj1);
  int nbin = h1->GetMaxSize();
  for(int j=0 ; j < nbin ; j++){
    double x;
    double y;
    double xerrhi = h1->GetErrorXhigh(j);
    double xerrlo = h1->GetErrorXlow(j);
    double yerrhi = h1->GetErrorYhigh(j);
    double yerrlo = h1->GetErrorYlow(j);
    int eff =  h1->GetPoint(j,x,y);
    cout << "[" << x-xerrlo  << "," << x+xerrhi << "] "  << " eff (" << ibin << ") = " << y << " (+" << yerrhi << " -" << yerrlo << ")" << endl;
    gr_Isomc->SetPoint(j,x,y);
    gr_Isomc->SetPointError(j,xerrlo,xerrhi,yerrlo,yerrhi);          
  }
  
  /////////////////////////////////////////////////////////
  cIso=new TCanvas("cIso","cIso",700,500);
  gPad->SetFillColor(0);

  gr_Iso->GetXaxis()->SetLimits(Xmin,Xmax);
  gr_Iso->SetMaximum(1.1);  
  gr_Iso->SetMinimum(hmin);
  gr_Iso->GetXaxis()->SetTitle(Form("%s",xtitle.Data()));
  gr_Iso->GetYaxis()->SetTitle(Form("%s",ytitle.Data()));

  gr_Iso->SetLineColor(4);
  gr_Iso->SetMarkerColor(4);
  gr_Iso->SetLineWidth(2);
  //gr_Iso->SetMarkerSize(1.1);
  gr_Iso->SetMarkerStyle(20);

  gr_Isomc->SetLineColor(2);
  gr_Isomc->SetMarkerColor(2);
  gr_Isomc->SetLineWidth(2);
  //gr_Isomc->SetMarkerSize(1.1);
  gr_Isomc->SetMarkerStyle(24);

  gr_Iso->Draw("APZ"); 
  gr_Isomc->Draw("PZsame"); 
  //TLegend *lIso= new TLegend(0.83,0.82,0.95,0.92);
  TLegend *lIso= new TLegend(0.5,0.2,0.75,0.3);
  lIso->SetBorderSize(0);
  lIso->SetTextFont(42);
  lIso->SetTextSize(0.04);
  lIso->SetLineColor(0);
  lIso->SetLineStyle(1);
  lIso->SetLineWidth(0.5);
  lIso->SetFillColor(0);
  lIso->SetFillStyle(1001);

  lIso->AddEntry(gr_Isomc,"  MC","PL");
  lIso->AddEntry(gr_Iso,"  Data","PL");

  lIso->Draw();



}
