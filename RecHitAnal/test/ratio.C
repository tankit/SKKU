#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <map>
#include <fstream.h>

void ratio()
{
  gROOT->LoadMacro("./tdrStyle.C");
  setTDRStyle();

  gStyle->SetPalette(1);
  
//  TFile* f_wRPC_RD = new TFile("Eff_data_medium_RPC.root");
//  TFile* f_woRPC_RD = new TFile("Eff_data_medium_woRPC.root");

//  TFile* f_wRPC_MC = new TFile("Eff_DY_medium_RPC.root");
//  TFile* f_woRPC_MC = new TFile("Eff_DY_medium_woRPC.root");

  TFile* f_wRPC_RD = new TFile("Eff_data_loose_RPC.root");
  TFile* f_woRPC_RD = new TFile("Eff_data_loose_woRPC.root");

  TFile* f_wRPC_MC = new TFile("Eff_DY_loose_RPC.root");
  TFile* f_woRPC_MC = new TFile("Eff_DY_loose_woRPC.root");

  TCanvas* c_wRPC_RD = (TCanvas*)f_wRPC_RD->Get("muonEffs/eta/fit_eff_plots/eta_PLOT");
  TCanvas* c_woRPC_RD = (TCanvas*)f_woRPC_RD->Get("muonEffs/eta/fit_eff_plots/eta_PLOT");

  TCanvas* c_wRPC_MC = (TCanvas*)f_wRPC_MC->Get("muonEffs/eta/fit_eff_plots/eta_PLOT");
  TCanvas* c_woRPC_MC = (TCanvas*)f_woRPC_MC->Get("muonEffs/eta/fit_eff_plots/eta_PLOT");

  RooHist* h_wRPC_RD = (RooHist*) c_wRPC_RD->FindObject("hxy_fit_eff");
  RooHist* h_woRPC_RD = (RooHist*) c_woRPC_RD->FindObject("hxy_fit_eff");

  RooHist* h_wRPC_MC = (RooHist*) c_wRPC_MC->FindObject("hxy_fit_eff");
  RooHist* h_woRPC_MC = (RooHist*) c_woRPC_MC->FindObject("hxy_fit_eff");
  
  if ( h_wRPC_RD->GetMaxSize() != h_woRPC_RD->GetMaxSize() )
  {
    cout << "Bin size is different" << endl;
    return;
  }
///////////////////////////////(Data)////////////////////////////////////////
///////////////////////////(with RPC)////////////////////////////////
  TGraphAsymmErrors* gr_wRPC_RD = new TGraphAsymmErrors();  

  const int nBin = h_wRPC_RD->GetMaxSize();
  for(int i=0; i < nBin; i++)
  {
    double x,y;
    double xerrhi = h_wRPC_RD->GetErrorXhigh(i);
    double xerrlo = h_wRPC_RD->GetErrorXlow(i);
    double yerrhi = h_wRPC_RD->GetErrorYhigh(i);
    double yerrlo = h_wRPC_RD->GetErrorYlow(i);
    int ibin =  h_wRPC_RD->GetPoint(i,x,y);
    gr_wRPC_RD->SetPoint(i,x,y);
    gr_wRPC_RD->SetPointError(i,xerrlo,xerrhi,yerrlo,yerrhi);

  }
//////////////////////////(without RPC)//////////////////////////////
  TGraphAsymmErrors* gr_woRPC_RD = new TGraphAsymmErrors();

  //const int nBin = h_woRPC_RD->GetMaxSize();
  for(int i=0; i < nBin; i++)
  {
    double x,y;
    double xerrhi = h_woRPC_RD->GetErrorXhigh(i);
    double xerrlo = h_woRPC_RD->GetErrorXlow(i);
    double yerrhi = h_woRPC_RD->GetErrorYhigh(i);
    double yerrlo = h_woRPC_RD->GetErrorYlow(i);
    int ibin =  h_woRPC_RD->GetPoint(i,x,y);
    gr_woRPC_RD->SetPoint(i,x,y);
    gr_woRPC_RD->SetPointError(i,xerrlo,xerrhi,yerrlo,yerrhi);

  }

  if ( h_wRPC_MC->GetMaxSize() != h_woRPC_MC->GetMaxSize() )
  {
    cout << "Bin size is different" << endl;
    return;
  }

///////////////////////////////(MC)////////////////////////////////////////
///////////////////////////(with RPC)////////////////////////////////
  TGraphAsymmErrors* gr_wRPC_MC = new TGraphAsymmErrors();

  const int nBin = h_wRPC_MC->GetMaxSize();
  for(int i=0; i < nBin; i++)
  {
    double x,y;
    double xerrhi = h_wRPC_MC->GetErrorXhigh(i);
    double xerrlo = h_wRPC_MC->GetErrorXlow(i);
    double yerrhi = h_wRPC_MC->GetErrorYhigh(i);
    double yerrlo = h_wRPC_MC->GetErrorYlow(i);
    int ibin =  h_wRPC_MC->GetPoint(i,x,y);
    gr_wRPC_MC->SetPoint(i,x,y);
    gr_wRPC_MC->SetPointError(i,xerrlo,xerrhi,yerrlo,yerrhi);

  }
//////////////////////////(without RPC)//////////////////////////////
  TGraphAsymmErrors* gr_woRPC_MC = new TGraphAsymmErrors();

  //const int nBin = h_woRPC_MC->GetMaxSize();
  for(int i=0; i < nBin; i++)
  {
    double x,y;
    double xerrhi = h_woRPC_MC->GetErrorXhigh(i);
    double xerrlo = h_woRPC_MC->GetErrorXlow(i);
    double yerrhi = h_woRPC_MC->GetErrorYhigh(i);
    double yerrlo = h_woRPC_MC->GetErrorYlow(i);
    int ibin =  h_woRPC_MC->GetPoint(i,x,y);
    gr_woRPC_MC->SetPoint(i,x,y);
    gr_woRPC_MC->SetPointError(i,xerrlo,xerrhi,yerrlo,yerrhi);

  }

/////////////////////////(Ratio_Data)////////////////////////////////
  TGraphAsymmErrors* gr_ratio_RD = new TGraphAsymmErrors();
  
  for(int i=0; i < nBin; i++)
  {
    double x_wRPC_RD, y_wRPC_RD;
    double x_woRPC_RD, y_woRPC_RD;
    double xerrlo_RD = gr_wRPC_RD->GetErrorXlow(i);
    double xerrhi_RD = gr_wRPC_RD->GetErrorXhigh(i);
    //double yerrlo_wRPC_RD = gr_wRPC_RD->GetErrorYlow(i);
    //double yerrhi_wRPC_RD = gr_wRPC_RD->GetErrorYhigh(i);
    //double yerrlo_woRPC_RD = gr_woRPC_RD->GetErrorYlow(i);
    //double yerrhi_woRPC_RD = gr_woRPC_RD->GetErrorYhigh(i);

    gr_wRPC_RD->GetPoint(i, x_wRPC_RD, y_wRPC_RD);
    gr_woRPC_RD->GetPoint(i, x_woRPC_RD, y_woRPC_RD);

    const double relErrlo_wRPC_RD = gr_wRPC_RD->GetErrorYlow(i)/y_wRPC_RD;
    const double relErrhi_wRPC_RD = gr_wRPC_RD->GetErrorYhigh(i)/y_wRPC_RD;
    const double relErrlo_woRPC_RD = gr_woRPC_RD->GetErrorYlow(i)/y_woRPC_RD;
    const double relErrhi_woRPC_RD = gr_woRPC_RD->GetErrorYhigh(i)/y_woRPC_RD;

    const double ratio_RD = (y_wRPC_RD == 0 ? 0 : y_woRPC_RD/y_wRPC_RD);
    //const double errlo = sqrt(relErrlo_wRPC_RD*relErrlo_wRPC_RD + relErrlo_woRPC_RD*relErrlo_woRPC_RD);
    const double errlo_RD = ratio_RD*TMath::Hypot(relErrlo_wRPC_RD, relErrlo_woRPC_RD);
    const double errhi_RD = ratio_RD*TMath::Hypot(relErrhi_wRPC_RD, relErrhi_woRPC_RD);

    gr_ratio_RD->SetPoint(i, x_wRPC_RD, ratio_RD);
    gr_ratio_RD->SetPointError(i, xerrlo_RD, xerrhi_RD, errlo_RD, errhi_RD);
  }
////////////////////////(Ratio_MC)////////////////////////////////
  TGraphAsymmErrors* gr_ratio_MC = new TGraphAsymmErrors();

  for(int i=0; i < nBin; i++)
  {
    double x_wRPC_MC, y_wRPC_MC;
    double x_woRPC_MC, y_woRPC_MC;
    double xerrlo_MC = gr_wRPC_MC->GetErrorXlow(i);
    double xerrhi_MC = gr_wRPC_MC->GetErrorXhigh(i);
    //double yerrlo_wRPC_MC = gr_wRPC_MC->GetErrorYlow(i);
    //double yerrhi_wRPC_MC = gr_wRPC_MC->GetErrorYhigh(i);
    //double yerrlo_woRPC_MC = gr_woRPC_MC->GetErrorYlow(i);
    //double yerrhi_woRPC_MC = gr_woRPC_MC->GetErrorYhigh(i);

    gr_wRPC_MC->GetPoint(i, x_wRPC_MC, y_wRPC_MC);
    gr_woRPC_MC->GetPoint(i, x_woRPC_MC, y_woRPC_MC);

    const double relErrlo_wRPC_MC = gr_wRPC_MC->GetErrorYlow(i)/y_wRPC_MC;
    const double relErrhi_wRPC_MC = gr_wRPC_MC->GetErrorYhigh(i)/y_wRPC_MC;
    const double relErrlo_woRPC_MC = gr_woRPC_MC->GetErrorYlow(i)/y_woRPC_MC;
    const double relErrhi_woRPC_MC = gr_woRPC_MC->GetErrorYhigh(i)/y_woRPC_MC;

    const double ratio_MC = (y_wRPC_MC == 0 ? 0 : y_woRPC_MC/y_wRPC_MC);
    //const double errlo_MC = sqrt(relErrlo_wRPC_MC*relErrlo_wRPC_MC + relErrlo_woRPC_MC*relErrlo_woRPC_MC);
    const double errlo_MC = ratio_MC*TMath::Hypot(relErrlo_wRPC_MC, relErrlo_woRPC_MC);
    const double errhi_MC = ratio_MC*TMath::Hypot(relErrhi_wRPC_MC, relErrhi_woRPC_MC);

    gr_ratio_MC->SetPoint(i, x_wRPC_MC, ratio_MC);
    gr_ratio_MC->SetPointError(i, xerrlo_MC, xerrhi_MC, errlo_MC, errhi_MC);
  }


  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  gr_ratio_RD->SetMaximum(1.1);
  gr_ratio_RD->SetMinimum(0.6);
  gr_ratio_RD->SetLineColor(kBlue);
  gr_ratio_RD->SetMarkerColor(kBlue);
  gr_ratio_RD->SetTitle(";Muon #eta;Ratio(withoutRPC/withRPC)");
  gr_ratio_RD->Draw("APZ");
//  gr_ratio_RD->GetXaxis()->SetRangeUser(-1.59,1.59);
  gr_ratio_RD->GetXaxis()->SetRangeUser(-1.8,1.8);
 
  gr_ratio_MC->SetMaximum(1.1);
  gr_ratio_MC->SetMinimum(0.6);
  gr_ratio_MC->SetLineColor(kRed);
  gr_ratio_MC->SetMarkerStyle(24);
  gr_ratio_MC->SetMarkerColor(kRed);
  gr_ratio_MC->Draw("PZ");
  gr_ratio_MC->GetXaxis()->SetRangeUser(-1.59,1.59);

  TLegend* legend = new TLegend(0.7,0.2,0.9,0.4);
  legend->AddEntry(gr_ratio_RD,"Data","LP");
  legend->AddEntry(gr_ratio_MC,"MC","LP");
  legend->SetFillColor(kWhite);
  legend->Draw();

  TLatex latex;
//  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  latex.DrawLatex(0.0,1.125,"Ratio_loose");
//  latex.DrawLatex(.3,.9,"K^{*0}");
//  latex.DrawLatex(.2,.8,longstring);

  

//  TFile* newfile = new TFile("newfile.root","recreate");
//  gr_ratio->Write();
//  newfile->Write();
//  newfile->Close();
}


