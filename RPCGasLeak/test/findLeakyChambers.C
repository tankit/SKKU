#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TMath.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

double temperature = 293;

bool hasValue(const std::vector<int>& arr, int value);

double func(double* x, double* p)
{
  return p[1]*x[0] + p[0];
}

void findLeakyChambers()
{
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(0.1);

  // Build Volume data map
  std::map<int, double> rpcVolumeMap;
  if ( true )
  {
    ifstream fin("../data/Volume.txt");
    double v;
    int chamberNumber;
    while ( !fin.eof() )
    {
      if ( !( fin >> v >> chamberNumber ) ) break;
      rpcVolumeMap[chamberNumber] = v;
    }
  }

  // Build leak data map
  const int nFile = 7;
  const int nColumn = 7;
  double minTime = 1e9, maxTime = -1e9;
  int minChamberNumber = 10000, maxChamberNumber = -1;
  std::map<int, TGraph*> rpcPressureMap;
  std::map<int, TGraph*> rpcGasLeakMap;

  for ( int i=0; i<nFile; ++i )
  {
    ifstream fin(Form("../data/DataSet%d.txt", i+1));
    int dummyI;
    fin >> dummyI; // dummy variable
    std::vector<int> chamberNumbers;
    for ( int j=0; j<nColumn; ++j )
    {
      int chamberNumber;
      fin >> chamberNumber;
      if ( minChamberNumber > chamberNumber ) minChamberNumber = chamberNumber;
      if ( maxChamberNumber < chamberNumber ) maxChamberNumber = chamberNumber;
      chamberNumbers.push_back(chamberNumber);
      if ( rpcPressureMap.find(chamberNumber) == rpcPressureMap.end() )
      {
        rpcPressureMap[chamberNumber] = new TGraph();
      }
      if ( rpcGasLeakMap.find(chamberNumber) == rpcGasLeakMap.end() )
      {
        rpcGasLeakMap[chamberNumber] = new TGraph();
      }
    }

    double time, mbar;

    while ( !fin.eof() )
    {
      if ( ! (fin >> time) ) break;
      if ( time < minTime ) minTime = time;
      if ( time > maxTime ) maxTime = time;
      for ( int j=0; j<nColumn; ++j )
      {
        fin >> mbar;
        const int chamberNumber = chamberNumbers[j];
        TGraph* grp = rpcPressureMap[chamberNumber];
        const int nPoint = grp->GetN();
        grp->SetPoint(nPoint, time, mbar);

        TGraph* grpLeak = rpcGasLeakMap[chamberNumber];
        grpLeak->SetPoint(nPoint, time, 273*mbar*rpcVolumeMap[chamberNumber]/temperature);
      }
    }
  }

  // Do the analysis
  const double threshold = 1.0;
  std::vector<int> leakyChambers;
  TF1* fLinear = new TF1("fLinear", func, minTime, maxTime, 2);
  //fLinear->SetParLimits(1, -1, 0);
  //const double fitMin1 = 900, fitMax1 = 900+1200;
  //const double fitMin2 = 900, fitMax2 = 900+1200;
  //const double fitMin3 = 900, fitMax3 = 900+1200;
  const double fitMin1 = 1500, fitMax1 = 1500+600;
  const double fitMin2 = 1200, fitMax2 = 1200+600;
  const double fitMin3 = 900, fitMax3 = 900+600;
  TGraphAsymmErrors* grpPressureDrop1 = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpPressureDrop2 = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpPressureDrop3 = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpLeakRate1 = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpLeakRate2 = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpLeakRate3 = new TGraphAsymmErrors();
  const double maxPressureDrop = 0.25;
  const double maxLeakRate = 100;
  TH1F* hPressureDrop1 = new TH1F("hPressureDrop1", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
  TH1F* hPressureDrop2 = new TH1F("hPressureDrop2", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
  TH1F* hPressureDrop3 = new TH1F("hPressureDrop3", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
  TH1F* hLeakRate1 = new TH1F("hLeakRate1", "Gas Leak Rate;Gas Leak Rate (mL/10 min)", 50, 0, maxLeakRate);
  TH1F* hLeakRate2 = new TH1F("hLeakRate2", "Gas Leak Rate;Gas Leak Rate (mL/10 min)", 50, 0, maxLeakRate);
  TH1F* hLeakRate3 = new TH1F("hLeakRate3", "Gas Leak Rate;Gas Leak Rate (mL/10 min)", 50, 0, maxLeakRate);

  for ( std::map<int, TGraph*>::const_iterator iter = rpcPressureMap.begin();
      iter != rpcPressureMap.end(); ++iter )
  {
    const int chamberNumber = iter->first;
    TGraph* grp = iter->second;
    TGraph* grpLeak = rpcGasLeakMap[chamberNumber];
    const int nData = grp->GetN();
    const double* xData = grp->GetX();
    const double* yData = grp->GetY();

    // Find leaky chambers
    if ( yData[nData-1] < threshold )
    {
      // Mark as leaky chamber
      leakyChambers.push_back(chamberNumber);
      cout << "Huge Leak chamber" << chamberNumber << endl;
      continue;
    }

    if ( true )
    {
      // Do fitting if it is not leaky chamber - Pressure Drop
      TFitResultPtr r1 = grp->Fit("fLinear", "QSEMW" , "", fitMin1, fitMax1);
      TFitResultPtr r2 = grp->Fit("fLinear", "+QSEMW", "", fitMin2, fitMax2);
      TFitResultPtr r3 = grp->Fit("fLinear", "+QSEMW", "", fitMin3, fitMax3);

      TList* functions = grp->GetListOfFunctions();
      TF1* f1 = (TF1*)functions->At(0);
      TF1* f2 = (TF1*)functions->At(1);
      TF1* f3 = (TF1*)functions->At(2);

      f1->SetLineColor(kRed);
      f2->SetLineColor(kGreen);
      f3->SetLineColor(kBlue);

      const int nGrpPressureDrop = grpPressureDrop1->GetN();
      //const double slope1 = r1->Value(1), errUp1 = r1->UpperError(1), errDn1 = r1->LowerError(1);
      //const double slope2 = r2->Value(1), errUp2 = r2->UpperError(1), errDn2 = r2->LowerError(1);
      //const double slope3 = r3->Value(1), errUp3 = r3->UpperError(1), errDn3 = r3->LowerError(1);
      const double slope1 = r1->Value(1), errUp1 = f1->GetParError(1), errDn1 = f1->GetParError(1);
      const double slope2 = r2->Value(1), errUp2 = f1->GetParError(1), errDn2 = f1->GetParError(1);
      const double slope3 = r3->Value(1), errUp3 = f1->GetParError(1), errDn3 = f1->GetParError(1);

      grpPressureDrop1->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope1);
      grpPressureDrop2->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope2);
      grpPressureDrop3->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope3);
      grpPressureDrop1->SetPointError(nGrpPressureDrop, 0.5, 0.5, 600*errDn1, 600*errUp1);
      grpPressureDrop2->SetPointError(nGrpPressureDrop, 0.5, 0.5, 600*errDn2, 600*errUp2);
      grpPressureDrop3->SetPointError(nGrpPressureDrop, 0.5, 0.5, 600*errDn3, 600*errUp3);

      hPressureDrop1->Fill(TMath::Min(maxPressureDrop-1e-9,-slope1*600));
      hPressureDrop2->Fill(TMath::Min(maxPressureDrop-1e-9,-slope2*600));
      hPressureDrop3->Fill(TMath::Min(maxPressureDrop-1e-9,-slope3*600));
    }

    if ( true )
    {
      // Do fitting if it is not leaky chamber - Leak Rate
      TFitResultPtr r1 = grpLeak->Fit("fLinear", "QSEMW" , "", fitMin1, fitMax1);
      TFitResultPtr r2 = grpLeak->Fit("fLinear", "+QSEMW", "", fitMin2, fitMax2);
      TFitResultPtr r3 = grpLeak->Fit("fLinear", "+QSEMW", "", fitMin3, fitMax3);

      TList* functions = grpLeak->GetListOfFunctions();
      TF1* f1 = (TF1*)functions->At(0);
      TF1* f2 = (TF1*)functions->At(1);
      TF1* f3 = (TF1*)functions->At(2);

      f1->SetLineColor(kRed);
      f2->SetLineColor(kGreen);
      f3->SetLineColor(kBlue);

      const int nGrpLeakRate = grpLeakRate1->GetN();
      //const double slope1 = r1->Value(1), errUp1 = r1->UpperError(1), errDn1 = r1->LowerError(1);
      //const double slope2 = r2->Value(1), errUp2 = r2->UpperError(1), errDn2 = r2->LowerError(1);
      //const double slope3 = r3->Value(1), errUp3 = r3->UpperError(1), errDn3 = r3->LowerError(1);
      const double slope1 = r1->Value(1), errUp1 = f1->GetParError(1), errDn1 = f1->GetParError(1);
      const double slope2 = r2->Value(1), errUp2 = f1->GetParError(1), errDn2 = f1->GetParError(1);
      const double slope3 = r3->Value(1), errUp3 = f1->GetParError(1), errDn3 = f1->GetParError(1);

      grpLeakRate1->SetPoint(nGrpLeakRate, chamberNumber, -600*slope1);
      grpLeakRate2->SetPoint(nGrpLeakRate, chamberNumber, -600*slope2);
      grpLeakRate3->SetPoint(nGrpLeakRate, chamberNumber, -600*slope3);
      grpLeakRate1->SetPointError(nGrpLeakRate, 0.5, 0.5, 600*errUp1, 600*errDn1);
      grpLeakRate2->SetPointError(nGrpLeakRate, 0.5, 0.5, 600*errUp2, 600*errDn2);
      grpLeakRate3->SetPointError(nGrpLeakRate, 0.5, 0.5, 600*errUp3, 600*errDn3);

      hLeakRate1->Fill(TMath::Min(maxLeakRate-1e-9,-slope1*600)*(-1.0));
      hLeakRate2->Fill(TMath::Min(maxLeakRate-1e-9,-slope2*600)*(-1.0));
      hLeakRate3->Fill(TMath::Min(maxLeakRate-1e-9,-slope3*600)*(-1.0));
    }

  }

  // Draw results
  TCanvas* cPressureVsTime = new TCanvas("cPressureVsTime", "Pressure vs Time", 500, 500);
  TH1F* hPressureVsTime = new TH1F("hPressureVsTime", "Pressure vs Time;Time (s);Pressure (mbar)", 100, minTime, maxTime);
  hPressureVsTime->SetMinimum(0);
  hPressureVsTime->SetMaximum(10);
  hPressureVsTime->Draw();
  for ( std::map<int, TGraph*>::const_iterator iter = rpcPressureMap.begin();
      iter != rpcPressureMap.end(); ++iter )
  {
    const int chamberNumber = iter->first;
    if ( hasValue(leakyChambers, chamberNumber) ) continue;
    TGraph* grp = iter->second;
    grp->Draw("p");
  }
/*
  TCanvas* cGasLeakVsTime = new TCanvas("cGasLeakVsTime", "GasLeak vs Time", 500, 500);
  TH1F* hGasLeakVsTime = new TH1F("hGasLeakVsTime", "GasLeak vs Time;Time (s);GasLeak (mbar)", 100, minTime, maxTime);
  hGasLeakVsTime->SetMinimum(0);
  hGasLeakVsTime->SetMaximum(150);
  hGasLeakVsTime->Draw();
  for ( std::map<int, TGraph*>::const_iterator iter = rpcGasLeakMap.begin();
      iter != rpcGasLeakMap.end(); ++iter )
  {
    const int chamberNumber = iter->first;
    if ( hasValue(leakyChambers, chamberNumber) ) continue;
    TGraph* grp = iter->second;
    grp->Draw("p");
  }
*/

  TCanvas* cPressureDropVsChamber = new TCanvas("cPressureDropVsChamber", "Pressure drop vs chamber", 500, 500);
  TH1F* hPressureDropVsChamber = new TH1F("hPressureDropVsChamber", "Pressure drop vs chamber;Chamber index;Pressure drop (mbar/10 min)", (maxChamberNumber-minChamberNumber+1), minChamberNumber, maxChamberNumber+1);
  hPressureDropVsChamber->SetMinimum(-0.1);//(-1e-3*600);
  hPressureDropVsChamber->SetMaximum(0.1);//(1e-4);
  hPressureDropVsChamber->Draw();
  grpPressureDrop1->SetLineColor(kRed);
  grpPressureDrop2->SetLineColor(kGreen);
  grpPressureDrop3->SetLineColor(kBlue);
  grpPressureDrop1->Draw("PZ");
  grpPressureDrop2->Draw("PZ");
  grpPressureDrop3->Draw("PZ");
  TLegend* legPressureDropVsChamber = new TLegend(0.6, 0.17, 0.9, 0.4);
  legPressureDropVsChamber->SetFillStyle(0);
  legPressureDropVsChamber->SetBorderSize(0);
  legPressureDropVsChamber->AddEntry(grpPressureDrop1, "1500-2100", "l");
  legPressureDropVsChamber->AddEntry(grpPressureDrop2, "1200-1800", "l");
  legPressureDropVsChamber->AddEntry(grpPressureDrop3, "900-1500" , "l");
  legPressureDropVsChamber->Draw();

  TCanvas* cLeakRateVsChamber = new TCanvas("cLeakRate", "Leak Rate vs Chamber", 500, 500);
  TH1F* hLeakRateVsChamber = new TH1F("hLeakRateVsChamber", "Leak Rate vs Chamber;Chamber index;Leak Rate (mL/10 min)", (maxChamberNumber-minChamberNumber+1), minChamberNumber, maxChamberNumber+1);
  hLeakRateVsChamber->SetMinimum(-1);
  hLeakRateVsChamber->SetMaximum(5);
  hLeakRateVsChamber->Draw();
  grpLeakRate1->SetLineColor(kRed);
  grpLeakRate2->SetLineColor(kGreen);
  grpLeakRate3->SetLineColor(kBlue);
  grpLeakRate1->Draw("PZ");
  grpLeakRate2->Draw("PZ");
  grpLeakRate3->Draw("PZ");
  TLegend* legLeakRateVsChamber = new TLegend(0.63, 0.6, 0.9, 0.9);
  legLeakRateVsChamber->SetFillStyle(0);
  legLeakRateVsChamber->SetBorderSize(0);
  legLeakRateVsChamber->AddEntry(grpLeakRate1, "1500-2100", "l");
  legLeakRateVsChamber->AddEntry(grpLeakRate2, "1200-1800", "l");
  legLeakRateVsChamber->AddEntry(grpLeakRate3, "900-1500", "l");
  legLeakRateVsChamber->Draw();

  TCanvas* cPressureDrop = new TCanvas("cPressureDrop", "PressureDrop", 500, 500);
  hPressureDrop1->SetLineColor(kRed);
  hPressureDrop2->SetLineColor(kGreen);
  hPressureDrop3->SetLineColor(kBlue);
  hPressureDrop1->SetMaximum(50);
  hPressureDrop1->Draw();
  hPressureDrop2->Draw("same");
  hPressureDrop3->Draw("same");
  TLegend* legPressureDrop = new TLegend(0.6, 0.6, 0.9, 0.9);
  legPressureDrop->SetFillStyle(0);
  legPressureDrop->SetBorderSize(0);
  legPressureDrop->AddEntry(hPressureDrop1, "1500-2100", "l");
  legPressureDrop->AddEntry(hPressureDrop2, "1200-1800", "l");
  legPressureDrop->AddEntry(hPressureDrop3, "900-1500" , "l");
  legPressureDrop->Draw();

  cPressureDropVsChamber->Print("cPressureDropVsChamber.png");
  cPressureDrop->Print("cPressureDrop.png");
  cPressureVsTime->Print("cPressureVsTime.png");
  cLeakRateVsChamber->Print("cLeakRateVsChamber.png");
}

bool hasValue(const std::vector<int>& arr, int value)
{
  for ( int i=0, n=arr.size(); i<n; ++i )
  {
    if ( arr[i] == value ) return true;
  }
  return false;
}

