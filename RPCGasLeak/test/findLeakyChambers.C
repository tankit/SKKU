#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

bool hasValue(const std::vector<int>& arr, int value);

void findLeakyChambers()
{
	gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(0.1);

	// Build Volume data map
	std::map<int, double> rpcVolumeMap;
	if ( true )
	{
		ifstream fin("../data/volume.txt");
		double v;
		for ( int i=0; i<=49; ++i )
		{
			fin >> v;
			rpcVolumeMap[i+1] = v;
		}
	}

	// Build leak data map
	const int nFile = 7;
	const int nColumn = 7;
	double minTime = 1e9, maxTime = -1e9;
	int minChamberNumber = 10000, maxChamberNumber = -1;
	std::map<int, TGraph*> rpcPressureMap;
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
			}
		}
	}

	// Do the analysis
	const double threshold = 1.0;
	std::vector<int> leakyChambers;
  TF1* fLinear = new TF1("fLinear", "[1]*x+[0]", minTime, maxTime);
  fLinear->SetParLimits(1, -1, 0);
	const double fitMin1 = 1500, fitMax1 = 1500+600;
	const double fitMin2 = 1200, fitMax2 = 1200+600;
	const double fitMin3 = 900, fitMax3 = 900+600;
	TGraphErrors* grpPressureDrop1 = new TGraphErrors();
	TGraphErrors* grpPressureDrop2 = new TGraphErrors();
	TGraphErrors* grpPressureDrop3 = new TGraphErrors();
	const double maxPressureDrop = 0.25;
	TH1F* hPressureDrop1 = new TH1F("hPressureDrop1", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
	TH1F* hPressureDrop2 = new TH1F("hPressureDrop2", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
	TH1F* hPressureDrop3 = new TH1F("hPressureDrop3", "Pressure Drop;Pressure drop (mbar/10 min)", 50, 0, maxPressureDrop);
	for ( std::map<int, TGraph*>::const_iterator iter = rpcPressureMap.begin();
		  iter != rpcPressureMap.end(); ++iter )
	{
		const int chamberNumber = iter->first;
		TGraph* grp = iter->second;
		const int nData = grp->GetN();
		const double* xData = grp->GetX();
		const double* yData = grp->GetY();

		// Find leaky chambers
		if ( yData[nData-1] < threshold )
		{
			// Mark as leaky chamber
			leakyChambers.push_back(chamberNumber);
			cout << "PressureDrop chamber" << chamberNumber << endl;
			continue;
		}

		// Do fitting if it is not leaky chamber
		grp->Fit("fLinear", "Q" , "", fitMin1, fitMax1);
		grp->Fit("fLinear", "+Q", "", fitMin2, fitMax2);
		grp->Fit("fLinear", "+Q", "", fitMin3, fitMax3);

		TList* functions = grp->GetListOfFunctions();
		TF1* f1 = (TF1*)functions->At(0);
		TF1* f2 = (TF1*)functions->At(1);
		TF1* f3 = (TF1*)functions->At(2);
	
		f1->SetLineColor(kRed);
		f2->SetLineColor(kGreen);
		f3->SetLineColor(kBlue);

		const int nGrpPressureDrop = grpPressureDrop1->GetN();
		double slope1 = f1->GetParameter(1), err1 = f1->GetParError(1);
		double slope2 = f2->GetParameter(1), err2 = f2->GetParError(1);
		double slope3 = f3->GetParameter(1), err3 = f3->GetParError(1);
		if ( slope1 > 0 ) slope1 = err1 = 0;
		if ( slope2 > 0 ) slope2 = err2 = 0;
		if ( slope3 > 0 ) slope3 = err3 = 0;
		grpPressureDrop1->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope1);
		grpPressureDrop2->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope2);
		grpPressureDrop3->SetPoint(nGrpPressureDrop, chamberNumber, 600*slope3);
		grpPressureDrop1->SetPointError(nGrpPressureDrop, 0, 600*err1);
		grpPressureDrop2->SetPointError(nGrpPressureDrop, 0, 600*err2);
		grpPressureDrop3->SetPointError(nGrpPressureDrop, 0, 600*err3);

		hPressureDrop1->Fill(TMath::Min(maxPressureDrop-1e-9,-slope1*600));
		hPressureDrop2->Fill(TMath::Min(maxPressureDrop-1e-9,-slope2*600));
		hPressureDrop3->Fill(TMath::Min(maxPressureDrop-1e-9,-slope3*600));
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

	TCanvas* cPressureDropVsChamber = new TCanvas("cPressureDropVsChamber", "Pressure drop vs chamber", 500, 500);
	TH1F* hPressureDropVsChamber = new TH1F("hPressureDropVsChamber", "Pressure drop vs chamber;Chamber index;Pressure drop (mbar/10 min)", (maxChamberNumber-minChamberNumber+1), minChamberNumber, maxChamberNumber+1);
	hPressureDropVsChamber->SetMinimum(-1e-3*600);
	hPressureDropVsChamber->SetMaximum(1e-4);
	hPressureDropVsChamber->Draw();
	grpPressureDrop1->SetLineColor(kRed);
	grpPressureDrop2->SetLineColor(kGreen);
	grpPressureDrop3->SetLineColor(kBlue);
	grpPressureDrop1->Draw("P");
	grpPressureDrop2->Draw("P");
	grpPressureDrop3->Draw("P");
	TLegend* legPressureDropVsChamber = new TLegend(0.6, 0.2, 0.9, 0.4);
	legPressureDropVsChamber->SetFillStyle(0);
	legPressureDropVsChamber->SetBorderSize(0);
	legPressureDropVsChamber->AddEntry(grpPressureDrop1, "1500-2100", "p");
	legPressureDropVsChamber->AddEntry(grpPressureDrop2, "1200-1800", "p");
	legPressureDropVsChamber->AddEntry(grpPressureDrop3, "900-1500" , "p");
	legPressureDropVsChamber->Draw();

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
}

bool hasValue(const std::vector<int>& arr, int value)
{
	for ( int i=0, n=arr.size(); i<n; ++i )
	{
		if ( arr[i] == value ) return true;
	}
	return false;
}

