// -*- C++ -*-
//
// Package:    TrackAssociator
// Class:      TestTrackAssociator
// 
/*

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Fri Apr 21 10:59:41 PDT 2006
// $Id: TestTrackAssociator.cc,v 1.5 2012/03/27 09:06:07 mskim Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

// calorimeter info
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include <boost/regex.hpp>

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "Utilities/Timing/interface/TimerStack.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQM/RPCMonitorDigi/interface/utils.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCChamber.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TrackingTools/TrackAssociator/interface/DetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TAMuonChamberMatch.h"

#include <TROOT.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <assert.h>
#include "TPaletteAxis.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TFolder.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#include "TNtuple.h"

using namespace edm;
using namespace std;
using namespace reco;

#define PI 3.141592654
#define TWOPI 6.283185308

static const int l = 7;
//static const int m = 10; //RE
static const int m = 13; //RE-, RE+

static const float maxEta = 1.8; //1.6; updated since 11/20/2011 //2.5; //updated 07/22/2011
static const float maxRes = 10; //cm
int nEvents = 0;

class TestTrackAssociator : public edm::EDAnalyzer {
 public:
   explicit TestTrackAssociator(const edm::ParameterSet&);
   virtual ~TestTrackAssociator(){
      TimingReport::current()->dump(std::cout);
   }
   
   virtual void analyze (const edm::Event&, const edm::EventSetup&);
   virtual void honeDivide(TH1* hdiv, const TH1* hnumer, const TH1* hdenom);
   virtual void htwoDivide(TH2* hdiv, const TH2* hnumer, const TH2* hdenom);
   virtual void beginJob();
   virtual void endJob();

 private:
   TrackDetectorAssociator trackAssociator_;
   TrackAssociatorParameters parameters_;

  RPCRecHitCollection* _ThePoints;
  RPCRecHitCollection* _TheRecHits;
  RPCRecHitCollection* _TheGlbHits;

  int maxNEvents_;
  int minPtTrk_;
  bool theStart;
  double rangestrips;
  string outputAscii_;

  // RPC efficiency map
  map<string,int> expected[3];
  map<string,int> rpcFound[3];
  map<string,int> rpcIndex;

  Int_t runNumber;
  Int_t eventNumber;
  Int_t nMuons;

  Double_t mean, rms;
  Double_t mean2, rms2;
  Int_t entries, nameId;
  Char_t name[40];

  // the Tree file
  //TTree *bxTree;
  //Int_t hitbx;
  //Double_t hitdx;
  Int_t denominator, numerator, detector;
  Double_t efficiency, error;

  TH2F *hresBX[13]; //0=All, 1~6:RB, 7~9:RE-, 10~12:RE+
  TH1F *hbigresBX[13], *hbigBXres[13]; 
  TH1F *hsmallresBX[13], *hsmallBXres[13];
 
  TH1F *hresRoll[2][2316], *hresMean[3], *hresRms[3], *hresEntries[3];
  TH2F *hresMeanRms[3], *hresMeanIndex[3], *hresRmsIndex[3];

  TH1F *hrechit[l], *hrechitRB[l], *hrechitRE[l], *hrechitRBRE[l];
  TH1F *hresidual[11], *hresidualRB[5], *hresidualRE[4], *hresidualRing[6], *hresidualRBlayer[7];
  TH1F *hres[11], *hresRB[5], *hresRE[4], *hresRing[6], *hresRBlayer[7];
  TH1F *hpull[11], *hpullRB[5], *hpullRE[4], *hpullRing[6], *hpullRBlayer[7];

  TH2F *htestRB[5], *htestRBlayer[7], *htestRE[4], *htestRing[6];

  TH1F *hbunchX[3], *hclusterSize[3];
  TH1F *heffRPCRecHit[3], *heffRPC[3], *herrRPC[3]; //RPC, RB, RE
  TH2F *heffErr[3];                                 //RPC, RB, RE
  TH1F *heffEta[3], *heffPhi[3]; //0=Crossed, 1=Matched, 2=Empty
  TH2F *heffEtaPhi[3];           //0=Crossed, 1=Matched, 2=Empty

  TH1F *hlocalXRB[3][5], *hlocalXLayer[3][7], *hlocalXRE[3][4], *htmpXRE[3][7];
  TH1F *hlocalXring[3][13], *hlocalXringTwo[3], *hlocalXringThree[3];

  TH1F *hpointNRB[2][5], *hpointNLayer[2][7], *hpointNRE[2][4], *htmpNRE[2][7];
  TH1F *hpointNring[2][13], *hpointNringTwo[2], *hpointNringThree[2];

  TH1F *hmuN;

  TH1F *hmupt[3][4], *hmueta[3][4], *hmuphi[3][4];
  TH1F *hpt[3][4][9], *heta[3][4][9], *hphi[3][4][9];

};

TestTrackAssociator::TestTrackAssociator(const edm::ParameterSet& iConfig)
{
   // TrackAssociator parameters
   edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
   parameters_.loadParameters( parameters );
   
   trackAssociator_.useDefaultPropagator();
   
   //cut
   minPtTrk_    = iConfig.getUntrackedParameter<int>("minPtTrk");
   maxNEvents_  = iConfig.getUntrackedParameter<int>("maxNEvents");
   rangestrips  = iConfig.getUntrackedParameter<double>("rangestrips",1.);
   outputAscii_ = iConfig.getParameter<string>("outputAscii");

   //produces<RPCRecHitCollection>("RPCRecHitAssociated");

   //now do what ever initialization is needed
   Service<TFileService> fs;

   // the rootple
   //bxTree = new TTree("bxTree","bxTree");

   hmuN      = fs->make<TH1F>("hmuN","N muons",10,0,10);

   for(int i=0; i<3; ++i) {
     int nbin = 100, xmax = 50; if(i==1) nbin = 100, xmax = 25; else if(i==2) nbin = 100, xmax = 5;
     hresMean[i]      = fs->make<TH1F>(Form("hresMean%d",i),"track-hit residuals Mean",nbin,-xmax,xmax);
     hresRms[i]       = fs->make<TH1F>(Form("hresRms%d",i),"track-hit residuals RMS",nbin,0,xmax*2);
     hresMeanRms[i]   = fs->make<TH2F>(Form("hresMeanRms%d",i),"track-hit residuals Mean vs Rms",nbin,0,xmax,nbin,-xmax,xmax);
     hresMeanIndex[i] = fs->make<TH2F>(Form("hresMeanIndex%d",i),"track-hit residuals Mean vs Index",2316,0,2316,nbin,-xmax,xmax);
     hresRmsIndex[i]  = fs->make<TH2F>(Form("hresRmsIndex%d",i),"track-hit residuals Rms vs Index",2316,0,2316,nbin,0,xmax);

     if(i==0) hresEntries[i] = fs->make<TH1F>(Form("hresEntries%d",i),"track-hit residuals Entries",100,0,100);
     else if(i==1) hresEntries[i] = fs->make<TH1F>(Form("hresEntries%d",i),"track-hit residuals Entries",100,0,200); 
     else if(i==2) hresEntries[i] = fs->make<TH1F>(Form("hresEntries%d",i),"track-hit residuals Entries",100,0,400);

     if(i<2) {
       int nbin2 = 100, xmax2 = 50; if(i==1) nbin2 = 200, xmax2 = 100; //if(i==1) nbin2 = xmax2 = 200;
       for(int j=0; j<2316; ++j) {
	 hresRoll[i][j] = fs->make<TH1F>(Form("hresRoll%d_%d",i,j),"track-hit residuals for roll",nbin2,-xmax2,xmax2);
       }
     }
   }

   int npt = 60, neta = 40, nphi = 50;

   for(int i=0; i<3; ++i) {
     for(int j=0; j<4; ++j) {

       hmupt[i][j]  = fs->make<TH1F>(Form("hmupt%d_%d",i,j) ,"mu pt ",npt,0,120);
       hmueta[i][j] = fs->make<TH1F>(Form("hmueta%d_%d",i,j),"mu eta",neta,-2.0,2.0);
       hmuphi[i][j] = fs->make<TH1F>(Form("hmuphi%d_%d",i,j),"mu phi",nphi,-PI,PI);

       if(i==2) {
         hmupt[i][j]->SetMarkerStyle(20); hmueta[i][j]->SetMarkerStyle(20); hmuphi[i][j]->SetMarkerStyle(20);
       }

       for(int k=0; k<9; ++k) {
         hpt[i][j][k]  = fs->make<TH1F>(Form("hpt%d_%d_%d",i,j,k) ,"mu pt ",npt,0,120);
         heta[i][j][k] = fs->make<TH1F>(Form("heta%d_%d_%d",i,j,k),"mu eta",neta,-2.0,2.0);
         hphi[i][j][k] = fs->make<TH1F>(Form("hphi%d_%d_%d",i,j,k),"mu phi",nphi,-PI,PI);

         if(i==2) {
           hpt[i][j][k]->SetMarkerStyle(20+k); heta[i][j][k]->SetMarkerStyle(20+k); hphi[i][j][k]->SetMarkerStyle(20+k);
           hpt[i][j][k]->SetMarkerColor(k+1); heta[i][j][k]->SetMarkerColor(k+1); hphi[i][j][k]->SetMarkerColor(k+1);
           hpt[i][j][k]->SetLineColor(k+1); heta[i][j][k]->SetLineColor(k+1); hphi[i][j][k]->SetLineColor(k+1);
         }

       }
     }
   }

   for(int j=0; j<13; ++j) {

     if(j<l) {
       hrechit[j]     = fs->make<TH1F>(Form("hrechit%d",j),  "nHits (Rpc)",12,-0.5,11.5);
       hrechitRE[j]   = fs->make<TH1F>(Form("hrechitRE%d",j),"nHits (RE) ",12,-0.5,11.5);
       hrechitRB[j]   = fs->make<TH1F>(Form("hrechitRB%d",j),"nHits (RB) ",12,-0.5,11.5);
       hrechitRBRE[j] = fs->make<TH1F>(Form("hrechitRBRE%d",j),"nHits (RB-RE) ",12,-0.5,11.5);
     }

     if(j<3) {
       hbunchX[j]       = fs->make<TH1F>(Form("hbunchX%d",j),"Bunch Crossing",11,-5.5,5.5);
       hclusterSize[j]  = fs->make<TH1F>(Form("hclusterSize%d",j),"Cluster Size",100,0,100);
       heffRPCRecHit[j] = fs->make<TH1F>(Form("heffRPCRecHit%d",j),"RecHits Efficiency",2,-0.5,1.5);
       //heffRPC[j]       = fs->make<TH1F>(Form("heffRPC%d",j),"Efficiency",110,0.,1.1);
       heffRPC[j]       = fs->make<TH1F>(Form("heffRPC%d",j),"Efficiency",55,0,1.1);
       herrRPC[j]       = fs->make<TH1F>(Form("herrRPC%d",j),"Efficiency error",50,0,0.5);
       heffErr[j]       = fs->make<TH2F>(Form("heffErr%d",j),"Efficiency vs. Error",50,0,0.5,55,0,1.1);
       heffRPC[j]->GetXaxis()->SetTitle("Efficiency"); heffRPC[j]->GetYaxis()->SetTitle("Number of rolls");
       herrRPC[j]->GetXaxis()->SetTitle("Efficiency error"); herrRPC[j]->GetYaxis()->SetTitle("Number of rolls");
       heffErr[j]->GetXaxis()->SetTitle("Error"); heffErr[j]->GetYaxis()->SetTitle("Efficiency");

       heffEta[j]       = fs->make<TH1F>(Form("heffEta%d",j),"RecHits Eta",neta,-2.0,2.0);
       heffPhi[j]       = fs->make<TH1F>(Form("heffPhi%d",j),"RecHits Phi",nphi,-PI,PI);
       heffEtaPhi[j]    = fs->make<TH2F>(Form("heffEtaPhi%d",j),"RecHits Eta-Phi",neta,-2.0,2.0,nphi,-PI,PI);

       if(j==2) {
	 heffEta[j]->SetMarkerStyle(20); heffEta[j]->SetMinimum(0); heffEta[j]->SetMaximum(1.1);
	 heffPhi[j]->SetMarkerStyle(20); heffPhi[j]->SetMinimum(0); heffPhi[j]->SetMaximum(1.1);
       }
     }

     int nbin = 100, xmax = 50;

     hresBX[j]    = fs->make<TH2F>(Form("hresBX%d",j),"Residual vs. BX",nbin,-xmax,xmax,11,-5.5,5.5);
     hbigresBX[j] = fs->make<TH1F>(Form("hbigresBX%d",j),"BX (|res|>30cm)",11,-5.5,5.5);
     hbigBXres[j] = fs->make<TH1F>(Form("hbigBXres%d",j),"residual (|BX|>0)",nbin,-xmax,xmax);
     hsmallresBX[j] = fs->make<TH1F>(Form("hsmallresBX%d",j),"BX (|res|<=30cm)",11,-5.5,5.5);
     hsmallBXres[j] = fs->make<TH1F>(Form("hsmallBXres%d",j),"residual (BX=0)",nbin,-xmax,xmax);

     if(j<11) hresidual[j]        = fs->make<TH1F>(Form("hresidual%d",j),"Residual",nbin,-xmax,xmax);
     if(j<5)  hresidualRB[j]      = fs->make<TH1F>(Form("hresidualRB%d",j),"Residual (RB)",nbin,-xmax,xmax);
     if(j<7)  hresidualRBlayer[j] = fs->make<TH1F>(Form("hresidualRBlayer%d",j),"Residual (RB)",nbin,-xmax,xmax);
     if(j<4)  hresidualRE[j]      = fs->make<TH1F>(Form("hresidualRE%d",j),"Residual (RE)",nbin,-xmax,xmax);
     if(j<6)  hresidualRing[j]    = fs->make<TH1F>(Form("hresidualRing%d",j),"Residual (Ring)",nbin,-xmax,xmax);

     nbin = 100, xmax = 25;
     if(j<11) hpull[j]        = fs->make<TH1F>(Form("hpull%d",j),"Pull",nbin,-xmax,xmax);
     if(j<5)  hpullRB[j]      = fs->make<TH1F>(Form("hpullRB%d",j),"Pull (RB)",nbin,-xmax,xmax);
     if(j<7)  hpullRBlayer[j] = fs->make<TH1F>(Form("hpullRBlayer%d",j),"Pull (RB)",nbin,-xmax,xmax);
     if(j<4)  hpullRE[j]      = fs->make<TH1F>(Form("hpullRE%d",j),"Pull (RE)",nbin,-xmax,xmax);
     if(j<6)  hpullRing[j]    = fs->make<TH1F>(Form("hpullRing%d",j),"Pull (Ring)",nbin,-xmax,xmax);

     int nbin2 = 200, xmax2 = 200;
     if(j<11) hres[j]        = fs->make<TH1F>(Form("hres%d",j),"Residual",nbin2,-xmax2,xmax2);
     if(j<5)  hresRB[j]      = fs->make<TH1F>(Form("hresRB%d",j),"Residual (RB)",nbin2,-xmax2,xmax2);
     if(j<7)  hresRBlayer[j] = fs->make<TH1F>(Form("hresRBlayer%d",j),"Residual (RB)",nbin2,-xmax2,xmax2);
     int nbin3 = 100, xmax3 = 100;
     if(j<4)  hresRE[j]      = fs->make<TH1F>(Form("hresRE%d",j),"Residual (RE)",nbin3,-xmax3,xmax3);
     if(j<6)  hresRing[j]    = fs->make<TH1F>(Form("hresRing%d",j),"Residual (Ring)",nbin3,-xmax3,xmax3);

     if(j<5) htestRB[j]      = fs->make<TH2F>(Form("htestRB%d",j),"dPhi (RB)",nphi,-PI,PI,nbin2,-PI,PI);
     if(j<7) htestRBlayer[j] = fs->make<TH2F>(Form("htestRBlayer%d",j),"dPhi (RB)",nphi,-PI,PI,nbin2,-PI,PI);
     if(j<4) htestRE[j]      = fs->make<TH2F>(Form("htestRE%d",j),"dPhi (RE)",nphi,-PI,PI,nbin,-PI,PI);
     if(j<6) htestRing[j]    = fs->make<TH2F>(Form("htestRing%d",j),"dPhi (Ring)",nphi,-PI,PI,nbin,-PI,PI);

     int nhit = 10;
     for(int i=0; i<3; ++i) {

       nbin = 130, xmax = 130;
       if(j<5) {
	 hlocalXRB[i][j]            = fs->make<TH1F>(Form("hlocalXRB%d_%d",i,j),"localX (RB)",nbin,-xmax,xmax);
	 if(i<2) hpointNRB[i][j]    = fs->make<TH1F>(Form("hpointNRB%d_%d",i,j),"pointN (RB)",nhit,0,nhit);
       }
       if(j<7) {
	 hlocalXLayer[i][j]         = fs->make<TH1F>(Form("hlocalXLayer%d_%d",i,j),"localX layer (RB)",nbin,-xmax,xmax);
	 if(i<2) hpointNLayer[i][j] = fs->make<TH1F>(Form("hpointNLayer%d_%d",i,j),"pointN layer (RB)",nhit,0,nhit);
       }

       nbin = 80, xmax = 80;
       if(j<4) {
	 hlocalXRE[i][j]         = fs->make<TH1F>(Form("hlocalXRE%d_%d",i,j),"localX (RE)",nbin,-xmax,xmax);
	 if(i<2) hpointNRE[i][j] = fs->make<TH1F>(Form("hpointNRE%d_%d",i,j),"pointN (RE)",nhit,0,nhit);
       }
       if(j<7) {
	 htmpXRE[i][j]           = fs->make<TH1F>(Form("htmpXRE%d_%d",i,j),"localX (RE)",nbin,-xmax,xmax);
	 if(i<2) htmpNRE[i][j]   = fs->make<TH1F>(Form("htmpNRE%d_%d",i,j),"pointN (RE)",nhit,0,nhit);
       }

       hlocalXring[i][j]         = fs->make<TH1F>(Form("hlocalXring%d_%d",i,j),"localX ring",nbin,-xmax,xmax);
       if(i<2) hpointNring[i][j] = fs->make<TH1F>(Form("hpointNring%d_%d",i,j),"pointN ring",nhit,0,nhit);
       if(j==0) {
         hlocalXringTwo[i]       = fs->make<TH1F>(Form("hlocalXringTwo%d",i),"localX ring",nbin,-xmax,xmax);
         hlocalXringThree[i]     = fs->make<TH1F>(Form("hlocalXringThree%d",i),"localX ring",nbin,-xmax,xmax);
	 if(i<2) {
	   hpointNringTwo[i]     = fs->make<TH1F>(Form("hpointNringTwo%d",i),"pointN ring",nhit,0,nhit);
	   hpointNringThree[i]   = fs->make<TH1F>(Form("hpointNringThree%d",i),"pointN ring",nhit,0,nhit);
	 }
       }

       if(i==2) {
	 hlocalXRB[i][j]->SetMarkerStyle(20); hlocalXRB[i][j]->SetMinimum(0); hlocalXRB[i][j]->SetMaximum(1.1);
	 hlocalXLayer[i][j]->SetMarkerStyle(20); hlocalXLayer[i][j]->SetMinimum(0); hlocalXLayer[i][j]->SetMaximum(1.1);
	 hlocalXRE[i][j]->SetMarkerStyle(20); hlocalXRE[i][j]->SetMinimum(0); hlocalXRE[i][j]->SetMaximum(1.1);
	 htmpXRE[i][j]->SetMarkerStyle(20); htmpXRE[i][j]->SetMinimum(0); htmpXRE[i][j]->SetMaximum(1.1);
	 hlocalXring[i][j]->SetMarkerStyle(20); hlocalXring[i][j]->SetMinimum(0); hlocalXring[i][j]->SetMaximum(1.1);
	 hlocalXringTwo[i]->SetMarkerStyle(20); hlocalXringTwo[i]->SetMinimum(0); hlocalXringTwo[i]->SetMaximum(1.1);
	 hlocalXringThree[i]->SetMarkerStyle(20); hlocalXringThree[i]->SetMinimum(0); hlocalXringThree[i]->SetMaximum(1.1);
       }

     }

   }

}

void TestTrackAssociator::honeDivide(TH1* hdiv, const TH1* hnumer, const TH1* hdenom) {

  for(int ibin=1; ibin<hnumer->GetNbinsX()+2; ++ibin) {
    float numer = hnumer->GetBinContent(ibin);
    float denom = hdenom->GetBinContent(ibin);
    float y = 0, yerr = 0;
    if(denom!=0) {
      y = numer/(double)denom;
      yerr = sqrt(y*(1-y)/(double)denom);
    }
    hdiv->SetBinContent(ibin,y); hdiv->SetBinError(ibin,yerr);
  }

}

void TestTrackAssociator::htwoDivide(TH2* hdiv, const TH2* hnumer, const TH2* hdenom) {

  for(int ibin=1; ibin<hnumer->GetNbinsX()+2; ++ibin) {
    for(int jbin=1; jbin<hnumer->GetNbinsY()+2; ++jbin) {

      float numer = hnumer->GetBinContent(ibin,jbin);
      float denom = hdenom->GetBinContent(ibin,jbin);
      float y = 0, yerr = 0;
      if(denom!=0) {
	y = numer/(double)denom;
	yerr = sqrt(y*(1-y)/(double)denom);
      }
      hdiv->SetBinContent(ibin,jbin,y); hdiv->SetBinError(ibin,jbin,yerr);
    }
  }

}

void TestTrackAssociator::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // select the event
   runNumber   = iEvent.id().run();
   eventNumber = iEvent.id().event();
   bool goodrun = runNumber==1 || (runNumber>163233 && runNumber<170249) || runNumber>170527 ;
   bool goodevt = (runNumber==1 || maxNEvents_==-1) ? true : nEvents < maxNEvents_ ;

   if( goodrun && goodevt ) {

     nEvents++;
     nMuons = 0;

     // get list of tracks and their vertices
     Handle<reco::MuonCollection> muons;
     iEvent.getByLabel("muons",muons);

     // get list of RPC hits
     Handle<RPCRecHitCollection > rpcRecHits;
     iEvent.getByLabel("rpcRecHits",rpcRecHits);
     
     // load RPC geometry
     ESHandle<RPCGeometry> rpcGeometry;
     iSetup.get<MuonGeometryRecord>().get(rpcGeometry);

     if (theStart) {  // define map with all rolls
       // load RPC geometry
       //ESHandle<RPCGeometry> rpcGeometry;
       //iSetup.get<MuonGeometryRecord>().get(rpcGeometry);
       vector< RPCRoll*> RPCrolls;
       RPCrolls = rpcGeometry->rolls() ;
       int nlines = 0;
       for(vector< RPCRoll*>::iterator RPCIt = RPCrolls.begin(); RPCIt != RPCrolls.end();RPCIt++){
         RPCDetId detId = (*RPCIt)->id();
         RPCGeomServ RPCname(detId);
	 rpcIndex[RPCname.name()] = nlines;
	 for(int i=0; i<3; ++i) {
	   expected[i][RPCname.name()] = 0; rpcFound[i][RPCname.name()] = 0;
	 }
         ///cout << rpcIndex[RPCname.name()] << " " << RPCname.name() << " " << detId.rawId() << endl;
	 nlines++;
       }
       theStart=false;
     }

     Handle<reco::BeamSpot> bsHandle;
     iEvent.getByLabel("offlineBeamSpot", bsHandle);
     const reco::TrackBase::Point & beamSpot = reco::TrackBase::Point(bsHandle->x0(), bsHandle->y0(), bsHandle->z0());
     /*
     Handle<vector<reco::Vertex> > pvHandle;
     iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
     //Handle<reco::VertexCollection> pvHandle;
     //iEvent.getByLabel(vertexSrc_, pvHandle);

     Handle<vector<reco::Vertex> > pvWithBsHandle;
     iEvent.getByLabel("offlinePrimaryVerticesWithBS", pvWithBsHandle);
     */

     //_ThePoints = new RPCRecHitCollection();
     edm::OwnVector<RPCRecHit> RPCPointVector;

     //_TheRecHits = new RPCRecHitCollection();
     edm::OwnVector<RPCRecHit> RPCHitVector;
     edm::OwnVector<RPCRecHit> GLBHitVector;

     vector<RPCRecHit*> rpcAllHitVector;
     //vector<RPCRecHit*> rpcGlbHitVector;
     //vector<TrackingRecHit*> rpcGlbHitVector;

     vector<DetId> rpcIdHit;
     vector<LocalPoint> rpcLocalPoints;

     //--http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_2_5/doc/html/d3/d45/HLTRPCTrigNoSyncFilter_8cc_source.html
     typedef RPCRecHitCollection::const_iterator RecHitIter;

     for (RecHitIter recHitIter = rpcRecHits->begin(); recHitIter != rpcRecHits->end(); ++recHitIter){
       if(!recHitIter->isValid()) continue;
       RPCDetId rollId = (RPCDetId)(*recHitIter).rpcId();
       LocalPoint recHitPos = recHitIter->localPosition();
       LocalError tmpErr    = recHitIter->localPositionError();
       const RPCRoll* rollassociated = rpcGeometry->roll(rollId);
       const BoundPlane & RPCSurface = rollassociated->surface(); 
       GlobalPoint RecHitInGlobal = RPCSurface.toGlobal(recHitPos);

       const int BX = recHitIter->BunchX(), cluSize = recHitIter->clusterSize();
       hbunchX[0]->Fill(BX); hclusterSize[0]->Fill(cluSize);
       //std::cout<<"\t \t We have an RPC Rec Hit! bx="<<BX<<" Roll="<<rollId<<" Global Position="<<RecHitInGlobal<<std::endl;

       //const int rawid = rollId.rawId(); // raw Id of the roll
       //rpcIdHit.push_back(rollId.rawId());
       rpcIdHit.push_back(rollId);
       rpcLocalPoints.push_back(recHitPos);

       RPCRecHit* rpcTmpHit = new RPCRecHit(rollId,BX,recHitPos,tmpErr);
       rpcAllHitVector.push_back(rpcTmpHit);
     }
     
     for(unsigned int i=0; i<rpcIdHit.size(); i++) {

       RPCDetId detId = rpcIdHit[i];
       //RPCGeomServ RPCname(detId);
       const RPCRoll * r1 = rpcGeometry->roll(detId);
       LocalPoint recHitLocal = rpcLocalPoints[i];
       GlobalPoint rpcGlobalPos = r1->surface().toGlobal(recHitLocal);

       const double globalX    = rpcGlobalPos.x();
       const double globalY    = rpcGlobalPos.y();
       const double globalZ    = rpcGlobalPos.z();
       const double globalPerp = rpcGlobalPos.perp();
       const double globalEta  = rpcGlobalPos.eta();
       const double globalPhi  = rpcGlobalPos.phi();
       /*
       LogVerbatim("TrackAssociator") << "Testing for roll "
				      << detId.rawId() << ", " //same to << roll->id() << ", "
				      << r1->id().region() << " ("
				      << r1->id().ring() << ", "
				      << r1->id().station() << ", "
				      << r1->id().sector() << ", "
				      << r1->id().layer() << ", "
				      << r1->id().subsector() << ", "
				      << r1->id().roll() << ") ";
       LogVerbatim("TrackAssociator") << "\t hit global point (x,y) and (z,perp,eta,phi): "
	  			      << globalX << ", "
	  			      << globalY << " and "
	  			      << globalZ << ", "
	 			      << globalPerp << ", "
	 			      << globalEta << ", "
	 			      << globalPhi ;
       LogVerbatim("TrackAssociator") << "\t hit local point (x,y): "
				      << recHitLocal.x() << ", "
				      << recHitLocal.y();
       */
     }
     
     //cut = cms.string(' pt > 20 && abs(eta)<2.4 && isGlobalMuon = 1 && isTrackerMuon = 1 && isolationR03().sumPt<3.0
     // && abs(innerTrack().dxy)<1.0
     // && innerTrack().numberOfValidHits()>10
     // && globalTrack().hitPattern().numberOfValidPixelHits()>0
     // && globalTrack().hitPattern().numberOfValidTrackerHits()>10
     // && globalTrack().normalizedChi2()<10.0
     // && globalTrack().hitPattern().numberOfValidMuonHits()>0 '),

     // loop 
     LogVerbatim("TrackAssociator") << "Number of muons found in the event: " << muons->size() ;
     for(reco::MuonCollection::const_iterator muon = muons->begin(); 
	 muon != muons->end(); ++muon){
       
       bool isTrackerMu = false, isStaMu = false, isGlobalMu = false;
       float ptMu = muon->pt(), etaMu = muon->eta(), phiMu = muon->phi();
       if (muon->isTrackerMuon() && abs(etaMu)<maxEta)    isTrackerMu = true;
       if (muon->isStandAloneMuon() && abs(etaMu)<maxEta) isStaMu     = true;
       if (muon->isGlobalMuon() && abs(etaMu)<maxEta)     isGlobalMu  = true;
       
       bool selMu = isGlobalMu && isTrackerMu && abs(etaMu)<maxEta && ptMu>=20;
       if (!selMu) {
         LogVerbatim("TrackAssociator") << "Skipped un-selected muon (Pt: " << muon->pt() << ", Eta: " << muon->eta() << ", Phi: " << muon->phi()
                                        << ", isGlobalMu: " << isGlobalMu << ", isTrackerMu: " << isTrackerMu << ", isStaMu: " << isStaMu << ")" ;
         continue;
       }

       //float d0 = -muon->innerTrack()->dxy();
       float dBs = -muon->innerTrack()->dxy(beamSpot); //or innerTrack()->dxy(bsHandle->position());
       /*
       float dBsError = (sqrt(pow(((*muons)[muidx]).innerTrack()->dxyError(),2)+pow(bsHandle->BeamWidthX(),2)+ pow(bsHandle->BeamWidthY(),2)));
       float dBsWidthX = bsHandle->BeamWidthX(), dBsWidthY = bsHandle->BeamWidthY();
       float dxy, dxyWithBs;
       math::XYZPoint vertexPosition(NAN, NAN, NAN);
       if (pvHandle.isValid()) {
         vertexPosition = (*pvHandle).size()>0 ? (*pvHandle)[0].position() : math::XYZPoint(NAN,NAN,NAN);
         dxy = ((*muons)[muidx]).innerTrack()->dxy(vertexPosition);
         //cout << " ===> vertexPosition    = " << vertexPosition << endl;
       }
       //--calculate dxy only in relation to the first reconstructed primary vertex
       std::vector<reco::Vertex>::const_iterator v = pvWithBsHandle->begin();
       if (v != pvWithBsHandle->end()) {
         const reco::Vertex& recoVertex = *v;
         dxyWithBs = muon->innerTrack()->dxy(recoVertex.position());
         //cout << " ===> vertexPositionTmp = " << recoVertex.position() << endl;
       }
       */

       float isoR03 = muon->isolationR03().sumPt;
       float normChi2 = muon->globalTrack()->normalizedChi2();
       int numValidHits = muon->innerTrack()->numberOfValidHits();
       int numValidMuonHits = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
       int numValidPixelHits = muon->globalTrack()->hitPattern().numberOfValidPixelHits();
       int numValidTrackerHits = muon->globalTrack()->hitPattern().numberOfValidTrackerHits();
       //cout << " ===> d0=" << d0 << " (|d0|=" << abs(d0) << "), isoR03=" << isoR03 << ", normChi2=" << normChi2 << ", nValidHits=" << numValidHits << ", nValidMuonHits=" << numValidMuonHits << endl;

       //bool selMu2 = abs(d0)<0.2 && isoR03<3.0 && normChi2<10.0 && numValidHits>10 && numValidMuonHits>0;
       //--https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatRecoVertex
       //--Very important update on May 17: remove d0 cut for test purpose because the simulated d0 (from 0,0,0) is much different from the real d0 in Data
       //if(Run==1) selMu2 = isoR03<3.0 && normChi2<10.0 && numValidHits>10 && numValidMuonHits>0;
       //--Don't forget to update selMuTest2 if selMu2 is updated !!!
       bool selMu2 = abs(dBs)<0.2 && isoR03<3.0 && normChi2<10.0 && numValidHits>10 && numValidMuonHits>0 && numValidPixelHits>0 && numValidTrackerHits>10;
 
       if (!selMu2) {
         LogVerbatim("TrackAssociator") << "Failed quality cuts (Pt: " << muon->pt() << ", Eta: " << muon->eta() << ")" ;
         continue;
       }

       TrackRef glbOfGlobalRef = muon->globalTrack();

       //if ( _ThePoints ) delete _ThePoints;
       //if ( _TheRecHits ) delete _TheRecHits;
       //if ( _TheGlbHits ) delete _TheGlbHits;

       _ThePoints  = new RPCRecHitCollection();
       _TheRecHits = new RPCRecHitCollection();
       _TheGlbHits = new RPCRecHitCollection();

       int hitsFromRpc = 0, hitsFromRB = 0, hitsFromRE = 0;
       for(trackingRecHit_iterator recHit = ((reco::Track)*glbOfGlobalRef).recHitsBegin(); recHit != ((reco::Track)*glbOfGlobalRef).recHitsEnd(); ++recHit){
	 if (!(*recHit)->isValid()) continue;
	 
	 DetId id = (*recHit)->geographicalId();
	 // hits from RPC
	 if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::RPC ) {
	   
	   const RPCRoll * r1;
	   r1 = rpcGeometry->roll(id);
	   RPCDetId rpcId = r1->id(); //Myrollid
	   LocalPoint recHitLocal = (*recHit)->localPosition();
	   LocalError tmpErr = (*recHit)->localPositionError();
	   GlobalPoint rpcGlobalPos= r1->surface().toGlobal(recHitLocal);
	   
	   const int region = rpcId.region();

           hitsFromRpc++; if(region==0) hitsFromRB++; else hitsFromRE++;
	   
	   //--/afs/cern.ch/user/m/mskim/testarea/CMSSW_3_8_4_patch2/src/RecoLocalMuon/RPCRecHit/src/RPCRecHitBaseAlgo.cc
	   //int bunchCrossing = 0; //cl->bx();
	   //int firstClustStrip = 0; //cl->firstStrip();
	   //int clusterSize = 0; //cl->clusterSize();
	   
	   RPCRecHit* rpcTmpHit = new RPCRecHit(rpcId,0,recHitLocal,tmpErr);
	   GLBHitVector.clear();
	   GLBHitVector.push_back(rpcTmpHit);
	   _TheGlbHits->put(rpcId,GLBHitVector.begin(),GLBHitVector.end());

	 }
       } //end loop on glbOfGlobalRef

       // skip low Pt tracks
       if (muon->pt() < 2) {
	 LogVerbatim("TrackAssociator") << "Skipped low Pt muon (Pt: " << muon->pt() << ")" ;
	 continue;
       }
      
       // skip tracks originated away from the IP
       if (fabs(muon->vertex().rho()) > 50) {
	 LogVerbatim("TrackAssociator") << "Skipped track with large impact parameter: " <<muon->vertex().rho();
	 continue;
       }
      
       LogVerbatim("TrackAssociator") << "\n-------------------------------------------------------\n Track (pt,eta,phi): " 
				      << muon->pt() << " , " << muon->eta() << " , " << muon->phi() ;

       //--For analysis     
       nMuons++;
       hrechit[0]->Fill(hitsFromRpc); hrechitRB[0]->Fill(hitsFromRB); hrechitRE[0]->Fill(hitsFromRE);
       if(abs(etaMu)<0.8) hrechitRB[1]->Fill(hitsFromRpc); else if(abs(etaMu)>1.2) hrechitRE[1]->Fill(hitsFromRpc); else hrechitRBRE[1]->Fill(hitsFromRpc);
       if(abs(etaMu)<0.8) hrechitRB[6]->Fill(hitsFromRB); else if(abs(etaMu)>1.2) hrechitRE[6]->Fill(hitsFromRE); else hrechitRBRE[6]->Fill(hitsFromRpc); //test

       int NrpcGlbHit = GLBHitVector.size(), NrpcAllHit = rpcAllHitVector.size();
       ///printf("Testing for event %10d: Pt=%6.2f, Eta=%5.2f with hitsFromRpc=%2d (cumulative: rpcGlbHit=%2d, rpcAllHit=%2d)\n",eventNumber,ptMu,etaMu,hitsFromRpc,NrpcGlbHit,NrpcAllHit);

       hmupt[0][0]->Fill(ptMu); hmueta[0][0]->Fill(etaMu); hmuphi[0][0]->Fill(phiMu);
       if(abs(etaMu)<0.8) { hmupt[0][1]->Fill(ptMu); hmuphi[0][1]->Fill(phiMu); }
       else if(abs(etaMu)>1.2) { hmupt[0][3]->Fill(ptMu); hmuphi[0][3]->Fill(phiMu); }
       else { hmupt[0][2]->Fill(ptMu); hmuphi[0][2]->Fill(phiMu); }

       for(int k=0; k<9; ++k) {
         if(k>hitsFromRpc) continue;
         hpt[0][0][k]->Fill(ptMu); heta[0][0][k]->Fill(etaMu); hphi[0][0][k]->Fill(phiMu);
         if(abs(etaMu)<0.8) { hpt[0][1][k]->Fill(ptMu); hphi[0][1][k]->Fill(phiMu); }
         else if(abs(etaMu)>1.2) { hpt[0][3][k]->Fill(ptMu); hphi[0][3][k]->Fill(phiMu); }
         else { hpt[0][2][k]->Fill(ptMu); hphi[0][2][k]->Fill(phiMu); }
       } 

       TrackDetMatchInfo info;
       if (muon->innerTrack().isAvailable()) 
	 info = trackAssociator_.associate(iEvent, iSetup, *muon->innerTrack(), parameters_);
       else {
	 if (!muon->outerTrack().isAvailable()) {
	   LogVerbatim("TrackAssociator") << "No refernced tracks are available, skim the muon";
	   continue;
	 }
	 info = trackAssociator_.associate(iEvent, iSetup, *muon->outerTrack(), parameters_);
       }
	   
       LogVerbatim("TrackAssociator") << "===========================================================================" ;
       LogVerbatim("TrackAssociator") << "ECAL RecHit energy: crossed, 3x3(max), 5x5(max), 3x3(direction), 5x5(direction), cone R0.5, generator";
       DetId centerId = info.findMaxDeposition(TrackDetMatchInfo::EcalRecHits);
       LogVerbatim("TrackAssociator") << "     " << 
	 info.crossedEnergy(TrackDetMatchInfo::EcalRecHits) << ", \t" <<
	 info.nXnEnergy(centerId, TrackDetMatchInfo::EcalRecHits, 1) << ", \t" <<
	 info.nXnEnergy(centerId, TrackDetMatchInfo::EcalRecHits, 2) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::EcalRecHits, 1) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::EcalRecHits, 2) << ", \t" <<
	 info.coneEnergy(0.5, TrackDetMatchInfo::EcalRecHits) << ", \t" <<
	 info.ecalTrueEnergy;
       LogVerbatim("TrackAssociator") << "ECAL trajectory point (z,Rho,eta,phi), max deposit DetId";
       LogVerbatim("TrackAssociator") << "     " <<
	 "(" << info.trkGlobPosAtEcal.z() << ", " << info.trkGlobPosAtEcal.Rho() << ", " <<
	 info.trkGlobPosAtEcal.eta() << ", " << info.trkGlobPosAtEcal.phi() << "), " << centerId.rawId();
       LogVerbatim("TrackAssociator") << "ECAL crossed DetIds with associated hits: (id, energy, z, perp, eta, phi)";
       for(std::vector<const EcalRecHit*>::const_iterator hit = info.crossedEcalRecHits.begin(); 
	   hit != info.crossedEcalRecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "ECAL crossed DetIds: (id, z, perp, eta, phi)";
       for(std::vector<DetId>::const_iterator id = info.crossedEcalIds.begin(); 
	   id != info.crossedEcalIds.end(); ++id)
	 {
	   GlobalPoint point = info.getPosition(*id);
	   LogVerbatim("TrackAssociator") << "\t" << id->rawId() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "ECAL associated DetIds: (id, energy, z, perp, eta, phi)";
       for(std::vector<const EcalRecHit*>::const_iterator hit = info.ecalRecHits.begin(); 
	   hit != info.ecalRecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "Preshower crossed DetIds: (id, z, perp, eta, phi)";
       for(std::vector<DetId>::const_iterator id = info.crossedPreshowerIds.begin(); 
	   id != info.crossedPreshowerIds.end(); ++id)
	 {
	   GlobalPoint point = info.getPosition(*id);
	   LogVerbatim("TrackAssociator") << "\t" << id->rawId() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }

       LogVerbatim("TrackAssociator") << "---------------------------------------------------------------------------" ;
       LogVerbatim("TrackAssociator") << "HCAL RecHit energy: crossed, 3x3(max), 5x5(max), 3x3(direction), 5x5(direction), cone R0.5, generator";
       centerId = info.findMaxDeposition(TrackDetMatchInfo::HcalRecHits);
       LogVerbatim("TrackAssociator") << "     " << 
	 info.crossedEnergy(TrackDetMatchInfo::HcalRecHits) << ", \t" <<
	info.nXnEnergy(centerId, TrackDetMatchInfo::HcalRecHits, 1) << ", \t" <<
	 info.nXnEnergy(centerId, TrackDetMatchInfo::HcalRecHits, 2) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::HcalRecHits, 1) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::HcalRecHits, 2) << ", \t" <<
	 info.coneEnergy(0.5, TrackDetMatchInfo::HcalRecHits) << ", \t" <<
	 info.hcalTrueEnergyCorrected;
       LogVerbatim("TrackAssociator") << "HCAL trajectory point (z,Rho,eta,phi), max deposit DetId";
       LogVerbatim("TrackAssociator") << "     " <<
	 "(" << info.trkGlobPosAtHcal.z() << ", " << info.trkGlobPosAtHcal.Rho() << ", " <<
	 info.trkGlobPosAtHcal.eta() << ", " << info.trkGlobPosAtHcal.phi() << "), " << centerId.rawId();
       LogVerbatim("TrackAssociator") << "HCAL crossed DetIds with hits:";
       for(std::vector<const HBHERecHit*>::const_iterator hit = info.crossedHcalRecHits.begin(); 
	   hit != info.crossedHcalRecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << (*hit)->id().depth() << ", \t" << 
	     point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "HCAL crossed DetIds: (id, z, perp, eta, phi)";
       for(std::vector<DetId>::const_iterator id = info.crossedHcalIds.begin(); 
	   id != info.crossedHcalIds.end(); ++id)
	 {
	   GlobalPoint point = info.getPosition(*id);
	   LogVerbatim("TrackAssociator") << "\t" << id->rawId() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "HCAL associated DetIds: id, (energy, z, perp, depth, eta, phi)";
       for(std::vector<const HBHERecHit*>::const_iterator hit = info.hcalRecHits.begin(); 
	   hit != info.hcalRecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << (*hit)->id().depth() << ", \t" << 
	     point.eta() << ", \t" << point.phi() << ")";
	 }
       
       LogVerbatim("TrackAssociator") << "---------------------------------------------------------------------------" ;
       LogVerbatim("TrackAssociator") << "HO RecHit energy: crossed, 3x3(max), 5x5(max), 3x3(direction), 5x5(direction), cone R0.5";
       centerId = info.findMaxDeposition(TrackDetMatchInfo::HORecHits);
       LogVerbatim("TrackAssociator") << "     " << 
	 info.crossedEnergy(TrackDetMatchInfo::HORecHits) << ", \t" <<
	 info.nXnEnergy(centerId, TrackDetMatchInfo::HORecHits, 1) << ", \t" <<
	 info.nXnEnergy(centerId, TrackDetMatchInfo::HORecHits, 2) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::HORecHits, 1) << ", \t" <<
	 info.nXnEnergy(TrackDetMatchInfo::HORecHits, 2) << ", \t" <<
	 info.coneEnergy(0.5, TrackDetMatchInfo::HORecHits);
       LogVerbatim("TrackAssociator") << "HO trajectory point (z,Rho,eta,phi), max deposit DetId";
       LogVerbatim("TrackAssociator") << "     " <<
	 "(" << info.trkGlobPosAtHO.z() << ", " << info.trkGlobPosAtHO.Rho() << ", " <<
	 info.trkGlobPosAtHO.eta() << ", " << info.trkGlobPosAtHO.phi() << "), " << centerId.rawId();
       LogVerbatim("TrackAssociator") << "HO crossed DetIds with hits:";
       for(std::vector<const HORecHit*>::const_iterator hit = info.crossedHORecHits.begin(); 
	   hit != info.crossedHORecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "HO crossed DetIds: (id, z, perp, eta, phi)";
       for(std::vector<DetId>::const_iterator id = info.crossedHOIds.begin(); 
	   id != info.crossedHOIds.end(); ++id)
	 {
	   GlobalPoint point = info.getPosition(*id);
	   LogVerbatim("TrackAssociator") << "\t" << id->rawId() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi() << ")";
	 }
       LogVerbatim("TrackAssociator") << "HO associated DetIds: (id, energy,position)";
       for(std::vector<const HORecHit*>::const_iterator hit = info.hoRecHits.begin(); 
	   hit != info.hoRecHits.end(); ++hit)
	 {
	   GlobalPoint point = info.getPosition((*hit)->detid());
	   LogVerbatim("TrackAssociator") << "\t" << (*hit)->detid().rawId() << ", " << (*hit)->energy() << 
	     " \t(" << point.z() << ", \t" << point.perp() << ", \t" << point.eta() << ", \t" << point.phi();
	   // << ") ## " << info.dumpGeometry((*hit)->detid());
	 }

       int hitFound[m]; for(int k=0; k<m; ++k) hitFound[k] = 0;

       if (parameters_.useMuon) {

         LogVerbatim("TrackAssociator") << "RPC matching details for event = " << eventNumber ;
         for(std::vector<TAMuonChamberMatch>::const_iterator rollassociated = info.chambers.begin();
             rollassociated!=info.chambers.end(); rollassociated++)
           {
             if (rollassociated->detector()!=3) continue;

	     RPCDetId rollId = rollassociated->id.rawId(); //same to rollassociated->id();
	     LocalPoint rpcLocalPos = rollassociated->tState.localPosition();
	     //LocalError rpcLocalErr = rollassociated->tState.localPositionError();

	     double localX = rollassociated->tState.localPosition().x();
	     double localDistX = rollassociated->localDistanceX ;
	     double localDistY = rollassociated->localDistanceY ;

	     const int region    = rollId.region();
	     const int ring      = rollId.ring();
	     const int sector    = rollId.sector();
	     const int station   = rollId.station();
	     const int sublayer  = rollId.layer();
	     const int subsector = rollId.subsector();
	     const int roll      = rollId.roll();
	     const int rawid     = rollId.rawId();
	     int layer = 0;
	     if (region==0) {
	       layer = station-1 + station*sublayer;
	       if ((station==2 && sublayer==2) || (station==4 && sublayer==1)) layer -= 1;
	     } else if(region==-1) layer = station + 6;
	     else layer = station + 9;

	     //const char* WHEEL = region==0 ? "wheel" : "ring ";
	     ///printf("-> Checking for event %10d: rollId=%d, region=%2d, wheel=%2d, station=%d, sector=%2d, layer=%d, subsector=%d, roll=%d, localX = %6.2f, edgeX = %6.2f, edgeY = %6.2f)\n",eventNumber,rawid,region,ring,station,sector,layer,subsector,roll,localX,localDistX,localDistY);

	     RPCRecHit RPCPoint(rollId,0,rpcLocalPos);
	     RPCPointVector.clear();
	     RPCPointVector.push_back(RPCPoint);
	     _ThePoints->put(rollId,RPCPointVector.begin(),RPCPointVector.end());

	     //--/afs/cern.ch/user/m/mskim/testarea/CMSSW_3_8_4_patch2/src/RecoLocalMuon/RPCRecHit/src/RPCRecHitAli.cc
             typedef std::pair<RPCRecHitCollection::const_iterator, RPCRecHitCollection::const_iterator> rangeRecHits;
             rangeRecHits recHitCollection =  rpcRecHits->get(rollassociated->id.rawId());
             //RPCRecHitCollection::const_iterator refHitIter;

	     RPCHitVector.clear();
	     for(RecHitIter refHitIter = recHitCollection.first; refHitIter != recHitCollection.second ; refHitIter++){
	       if(!refHitIter->isValid()) continue;

               hitFound[layer]++;
 
	       LocalPoint recHitPos = refHitIter->localPosition();
	       LocalError tmpErr    = refHitIter->localPositionError();

	       const int BX = refHitIter->BunchX(), cluSize = refHitIter->clusterSize();
	       hbunchX[1]->Fill(BX); hclusterSize[1]->Fill(cluSize);

	       RPCRecHit* rpcTmpHit = new RPCRecHit(rollId,BX,recHitPos,tmpErr);
	       RPCHitVector.push_back(rpcTmpHit);
	     }   
	     _TheRecHits->put(rollId,RPCHitVector.begin(),RPCHitVector.end());
	   }

	 


         /*
	 LogVerbatim("TrackAssociator") << "Muon detector matching details: " ;
	 for(std::vector<TAMuonChamberMatch>::const_iterator chamber = info.chambers.begin();
	     chamber!=info.chambers.end(); chamber++)
	   {
	     LogVerbatim("TrackAssociator") << chamber->info() << "\n\t(DetId, station, edgeX, edgeY): "
					    << chamber->id.rawId() << ", "
					    << chamber->station() << ", "
					    << chamber->localDistanceX << ", "
					    << chamber->localDistanceY << ", ";
	     LogVerbatim("TrackAssociator") << "\t trajectory global point (z,perp,eta,phi): "
					    << chamber->tState.globalPosition().z() << ", "
					    << chamber->tState.globalPosition().perp() << ", "
					    << chamber->tState.globalPosition().eta() << ", "
					    << chamber->tState.globalPosition().phi() ;
	     LogVerbatim("TrackAssociator") << "\t trajectory local point (x,y): "
					    << chamber->tState.localPosition().x() << ", "
					    << chamber->tState.localPosition().y();
	     
	     for(std::vector<TAMuonSegmentMatch>::const_iterator segment=chamber->segments.begin(); 
		 segment!=chamber->segments.end(); segment++)
	       {
		 LogVerbatim("TrackAssociator") << "\t segment position (z,Rho,eta,phi,DetId): " 
						<< segment->segmentGlobalPosition.z() << ", "
						<< segment->segmentGlobalPosition.Rho() << ", "
						<< segment->segmentGlobalPosition.eta() << ", "
						<< segment->segmentGlobalPosition.phi() << ", "
						<< chamber->id.rawId();
		 LogVerbatim("TrackAssociator") << "\t segment local position (x,y): "
						<< segment->segmentLocalPosition.x() << ", "
						<< segment->segmentLocalPosition.y();
	       }
	   }
           */
       } //MuonSystem


       // Start matching RefHits to RecHits
       typedef std::map<RecHitIter, RecHitIter> RecToRecHitMap;
       RecToRecHitMap refToRecHitMap;

       int nFound[m], nCrossed[m], nMatched[m], minDxRing[m];
       double minDx[m], minDphi[m];
       for(int k=0; k<m; ++k) { nFound[k] = nCrossed[k] = nMatched[k] = 0; minDxRing[k] = 0; minDx[k] = minDphi[k] = -999; }

       RPCRecHitCollection::const_iterator recHitIterTmp[m];
       RPCRecHitCollection::const_iterator refHitIterTmp[m];
       //std::vector<RPCRecHit>::const_iterator muonRecHitIter;

       for (RecHitIter recHitIter = _ThePoints->begin(); recHitIter != _ThePoints->end(); ++recHitIter) {
	 //if (!recHitIter->isValid()) continue;

	 const RPCDetId recDetId = static_cast<const RPCDetId>(recHitIter->rpcId());
	 const RPCRoll* recRoll = dynamic_cast<const RPCRoll*>(rpcGeometry->roll(recDetId));
	 if ( !recRoll ) continue;

         RPCGeomServ RPCname(recDetId);
	 const int region    = recDetId.region();
	 const int ring      = recDetId.ring();
	 const int sector    = recDetId.sector();
	 const int station   = recDetId.station();
	 const int sublayer  = recDetId.layer();
	 const int subsector = recDetId.subsector();
	 const int roll      = recDetId.roll();
	 const int rawid     = recDetId.rawId();
	 int layer = 0;
	 if (region==0) {
	   layer = station-1 + station*sublayer;
	   if ((station==2 && sublayer==2) || (station==4 && sublayer==1)) layer -= 1;
	 } else if(region==-1) layer = station + 6;
	 else layer = station + 9;
	 
         int ringTmp = 0;
         if(region==-1) {
           if(ring==2) {
             if(station==1) ringTmp = 1;       //RE ring -2
             else if(station==2) ringTmp = 2;  //RE ring -2
             else if(station==3) ringTmp = 3;  //RE ring -2
           } else if(ring==3) {
             if(station==1) ringTmp = 7;       //RE ring -3
             else if(station==2) ringTmp = 8;  //RE ring -3
             else if(station==3) ringTmp = 9;  //RE ring -3
           }
         } else if(region==1) {
           if(ring==2) {
             if(station==1) ringTmp = 4;       //RE ring +2
             else if(station==2) ringTmp = 5;  //RE ring +2
             else if(station==3) ringTmp = 6;  //RE ring +2
           } else if(ring==3) {
             if(station==1) ringTmp = 10;      //RE ring +3
             else if(station==2) ringTmp = 11; //RE ring +3
             else if(station==3) ringTmp = 12; //RE ring +3
           }
         }

	 LocalPoint recHitLocal   = recHitIter->localPosition();
	 GlobalPoint rpcGlobalPos = recRoll->surface().toGlobal(recHitLocal);
	 //const double globalZ     = rpcGlobalPos.z(), globalPerp = rpcGlobalPos.perp();
	 const double globalEta   = rpcGlobalPos.eta(), globalPhi = rpcGlobalPos.phi(); 
	 const double localX      = recHitLocal.x();

	 nCrossed[layer]++;
	 if(nCrossed[layer]==1) {
	   heffRPCRecHit[0]->Fill(0); if(region==0) heffRPCRecHit[1]->Fill(0); else heffRPCRecHit[2]->Fill(0);
	 }

	 //const char* WHEEL = region==0 ? "wheel" : "ring ";
	 ///printf("--> Checking for event %10d: rollId=%d, region=%2d, wheel=%2d, station=%d, sector=%2d, layer=%d, subsector=%d, roll=%d, localX = %6.2f)\n",eventNumber,rawid,region,ring,station,sector,layer,subsector,roll,localX);
	 ///cout << " | Crossed : " << RPCname.name() << " with nCrossed = " << nCrossed[layer] << " for event " << eventNumber << endl;

	 for (RecHitIter refHitIter = _TheRecHits->begin(); refHitIter != _TheRecHits->end(); ++refHitIter) {
	 //for (RecHitIter refHitIter = _TheGlbHits->begin(); refHitIter != _TheGlbHits->end(); ++refHitIter) {

	   //if (!refHitIter->isValid()) continue;

	   const RPCDetId refDetId = static_cast<const RPCDetId>(refHitIter->rpcId());
	   const RPCRoll* refRoll = dynamic_cast<const RPCRoll*>(rpcGeometry->roll(refDetId));
	   if ( !refRoll ) continue;

	   RPCGeomServ hitRPCname(refDetId);
	   const int hitRegion    = refDetId.region();
	   const int hitRing      = refDetId.ring();
	   const int hitSector    = refDetId.sector();
	   const int hitStation   = refDetId.station();
	   const int hitSublayer  = refDetId.layer();
	   const int hitSubsector = refDetId.subsector();
	   const int hitRoll      = refDetId.roll();
	   const int hitRawid     = refDetId.rawId();
	   int hitLayer = 0;
	   if (hitRegion==0) {
	     hitLayer = hitStation-1 + hitStation*hitSublayer;
	     if ((hitStation==2 && hitSublayer==2) || (hitStation==4 && hitSublayer==1)) hitLayer -= 1;
	   } else if(region==-1) hitLayer = hitStation + 6;
	   else hitLayer = hitStation + 9;

           //if ( hitLayer == layer && hitSector == sector && hitSubsector == subsector && hitRoll == roll) nFound[layer]++;

	   //int deltaSector = fabs(hitSector - sector);
	   //if (deltaSector>6) deltaSector = 12 - deltaSector; 
	   //if ( (hitLayer != layer) || deltaSector>1 ) continue;

	   //if ( (hitLayer != layer) ) continue;                                                   //test1
	   //if ( (hitLayer != layer) || (region != 0 && hitSubsector != subsector) ) continue;     //test2
	   //if ( (hitLayer != layer) || (hitSubsector != subsector) ) continue;                    //test3
	   if ( (hitLayer != layer) || (hitSubsector != subsector) || (hitRoll != roll) ) continue; //testing ok
	   //---double check
	   if ( region!=0 && hitRing != ring ) {
	     ///cout << " Double check for different ring " << ring << " vs. " << hitRing << " in region = " << region << endl;
	     continue;
	   }

	   LocalPoint hitLocal       = refHitIter->localPosition();
	   GlobalPoint hitGlobalPos  = refRoll->surface().toGlobal(hitLocal);
	   //const double hitGlobalZ   = hitGlobalPos.z(), hitGlobalPerp = hitGlobalPos.perp(), hitGlobalEta = hitGlobalPos.eta();
	   const double hitGlobalPhi = hitGlobalPos.phi();
	   const double hitLocalX    = hitLocal.x();

	   const double newDx   = localX - hitLocalX;
	   const double newDphi = globalPhi - hitGlobalPhi;

	   //--http://cmslxr.fnal.gov/lxr/source/HLTrigger/special/src/HLTRPCFilter.cc#080
	   const int BX = refHitIter->BunchX(), cluSize = refHitIter->clusterSize();
	   hbunchX[2]->Fill(BX); hclusterSize[2]->Fill(cluSize); //here, clusterSize is always 99 while bunchX is OK
	   hresBX[0]->Fill(newDx,BX); hresBX[layer]->Fill(newDx,BX);
	   if(fabs(newDx) > 30) { hbigresBX[0]->Fill(BX); hbigresBX[layer]->Fill(BX); }
           else { hsmallresBX[0]->Fill(BX); hsmallresBX[layer]->Fill(BX); } 
	   if(fabs(BX) > 0) { hbigBXres[0]->Fill(newDx); hbigBXres[layer]->Fill(newDx); }
           else if(BX==0) { hsmallBXres[0]->Fill(newDx); hsmallBXres[layer]->Fill(newDx); }

	   string nameString = RPCname.name();
	   sscanf(nameString.c_str(),"%s",name);
	   
	   //TString nameChamber = nameString;
	   //if(nameChamber.Contains("RB")) detector = 1;
	   //else if(nameChamber.Contains("RE")) detector = 2;

	   //hitbx = BX; hitdx = newDx;

	   //if(fabs(newDx)>=(rangestrips+cluSize*0.5)*3){ //3 is a typical strip width for RPCs
	   //  cout <<"RPC passed, RecHits but far away, (rangestrips+cluSize*0.5)*3 = "<< (rangestrips+cluSize*0.5)*3 << ") " << refDetId << endl;
	   //}

	   minDxRing[layer] = fabs(newDx) < fabs(minDx[layer]) ? ring  : minDxRing[layer];
	   minDx[layer]     = fabs(newDx) < fabs(minDx[layer]) ? newDx : minDx[layer];
	   minDphi[layer]   = fabs(newDphi) < fabs(minDphi[layer]) ? newDphi : minDphi[layer];
	   
	   //const char* hitWHEEL = hitRegion==0 ? "wheel" : "ring ";
	   ///printf("%22s rpcRecHits: rollId=%d, region=%2d, wheel=%2d, station=%d, sector=%2d, layer=%d, subsector=%d, roll=%d, localX = %6.2f)\n","",hitRawid,hitRegion,hitRing,hitStation,hitSector,hitLayer,hitSubsector,hitRoll,hitLocalX);
	   ///printf("==> Residual for event %10d: deltaX=%6.2f (localX = %6.2f, hitLocalX = %6.2f), dPhi=%f (globalPhi = %f, hitGlobalPhi = %f)\n",eventNumber,newDx,localX,hitLocalX,newDphi,globalPhi,hitGlobalPhi);
	   ///printf("===> RPCname = %s: dx = %6.2f, bx = %d\n",name,newDx,BX);

	   // Associate RefHit to RecHit
	   RecToRecHitMap::const_iterator prevRefToReco = refToRecHitMap.find(recHitIter);

	   if ( prevRefToReco == refToRecHitMap.end() ) {
	     nMatched[layer]++;
	     ///cout << " | Matched : " << hitRPCname.name() << " with nMatched = " << nMatched[layer] << " in region and layer (" << region << ", " << layer << ") for ptMu and etaMu (" << ptMu << ", " << etaMu << ")" << endl;
	     if ( nMatched[layer] == 1 ) {
	       refToRecHitMap.insert(std::make_pair(recHitIter, refHitIter));
	       recHitIterTmp[layer] = recHitIter;
	       refHitIterTmp[layer] = refHitIter;
	       //rpcFound[0][RPCname.name()]++; if(region==0) rpcFound[1][RPCname.name()]++; else rpcFound[2][RPCname.name()]++;
	       ///cout << " | Inserted: " << RPCname.name() << " " << hitRPCname.name() << endl;
	     }
	   } else {
	     const double oldDx = (prevRefToReco->first->localPosition().x() - prevRefToReco->second->localPosition().x());
	     
	     //cout << " Printing newDx = " << newDx << " " << " oldDx = " << oldDx << " |" << prevRefToReco->first->localPosition().x() << "-" << prevRefToReco->second->localPosition().x() << "|" << endl;
	     
	     if ( fabs(newDx) < fabs(oldDx) ) {
	       ///cout << "   New Vs old: newDx = " << fabs(newDx) << " " << " oldDx = " << fabs(oldDx) << endl;
	       ///cout << " | Updated : " << RPCname.name() << " " << hitRPCname.name() << endl;
	       refToRecHitMap[recHitIter] = refHitIter;
	       refHitIterTmp[layer] = refHitIter;
	     }
	   }
	   
	   if ( nMatched[layer] > 1 ) {
	     //RecToRecHitMap::const_iterator tempRefToReco = refToRecHitMap.find(recHitIterTmp[layer]);
	     //const double oldDx = (tempRefToReco->first->localPosition().x() - tempRefToReco->second->localPosition().x());
	     const double oldDx = (recHitIterTmp[layer]->localPosition().x() - refHitIterTmp[layer]->localPosition().x());
	     ///cout << " Double counting: newDx = " << fabs(newDx) << " " << " oldDx = " << fabs(oldDx) << endl;

	     if ( fabs(newDx) < fabs(oldDx) ) {
	       //refToRecHitMap.erase(tempRefToReco); //compile error
	       refToRecHitMap.erase(recHitIterTmp[layer]); //delete &refHitIter; //running error
	       ///cout << " | Erased : " << endl;
	       refToRecHitMap.insert(std::make_pair(recHitIter, refHitIter));
	       ///cout << " | Re-inserted : " << RPCname.name() << " " << hitRPCname.name() << endl;
	       recHitIterTmp[layer] = recHitIter;
	       refHitIterTmp[layer] = refHitIter;

	     }

	   }

	 }

         if (hitFound[layer]==0 && nCrossed[layer]>=1) {
	   
	   ///cout << " | Not inserted: " << RPCname.name() << " with nCrossed = " << nCrossed[layer] << " in layer = " << layer << endl;

	   if (nCrossed[layer]==1) {
	     
	     expected[0][RPCname.name()]++; if(region==0) expected[1][RPCname.name()]++; else expected[2][RPCname.name()]++;
	     heffEta[0]->Fill(globalEta); heffPhi[0]->Fill(globalPhi); heffEtaPhi[0]->Fill(globalEta,globalPhi);
	     
	     if (region==0) {
	       hlocalXRB[0][0]->Fill(localX); hlocalXRB[0][station]->Fill(localX); hlocalXLayer[0][0]->Fill(localX); hlocalXLayer[0][layer]->Fill(localX);
	     } else {
	       int diskTmp = layer - 6;
	       hlocalXRE[0][0]->Fill(localX); hlocalXRE[0][station]->Fill(localX); htmpXRE[0][0]->Fill(localX); htmpXRE[0][diskTmp]->Fill(localX);
	       hlocalXring[0][0]->Fill(localX); hlocalXring[0][ringTmp]->Fill(localX);
	       if(ring==2) hlocalXringTwo[0]->Fill(localX); else if(ring==3) hlocalXringThree[0]->Fill(localX);
	     }
	   }
	   
	 } else if(hitFound[layer]!=0) {
	   if (region==0) {
	     hpointNRB[1][0]->Fill(nCrossed[layer]); hpointNRB[1][station]->Fill(nCrossed[layer]); hpointNLayer[1][0]->Fill(nCrossed[layer]); hpointNLayer[1][layer]->Fill(nCrossed[layer]);
	   } else {
	     int diskTmp = layer - 6;
	     hpointNRE[1][0]->Fill(nCrossed[layer]); hpointNRE[1][station]->Fill(nCrossed[layer]); htmpNRE[1][diskTmp]->Fill(nCrossed[layer]);
	     hpointNring[1][0]->Fill(nCrossed[layer]); hpointNring[1][ringTmp]->Fill(nCrossed[layer]);
	     if(ring==2) hpointNringTwo[1]->Fill(nCrossed[layer]); else if(ring==3) hpointNringThree[1]->Fill(nCrossed[layer]);
	   }
	 }

       }

       // Now we have refHit-recHit mapping
       // So we can fill up relavant histograms
       cout << " " << endl << " Summary for run/event " << runNumber << " / " << eventNumber << ": nEvents = " << nEvents << endl;
       
       for (int k=1; k<m; ++k) {
	 if (hitFound[k]==0 && nCrossed[k]>=1) {
	   
	   int layer = k, station = 0, diskTmp = 0;

	   if (k<7) {
	     if(k<=2) station = 1; else if(k<=4) station = 2; else if(k==5) station = 3; else station = 4;
	     hpointNRB[0][0]->Fill(nCrossed[k]); hpointNRB[0][station]->Fill(nCrossed[k]); hpointNLayer[0][0]->Fill(nCrossed[k]); hpointNLayer[0][k]->Fill(nCrossed[k]); 
	   } else {
	     diskTmp = k - 6;
	     station = diskTmp; if(k>9) station -= 3;
	     hpointNRE[0][0]->Fill(nCrossed[k]); hpointNRE[0][station]->Fill(nCrossed[k]); htmpNRE[0][diskTmp]->Fill(nCrossed[k]);
	     //hpointNring[0][0]->Fill(nCrossed[k]); hpointNring[0][ringTmp]->Fill(nCrossed[k]);
	     //if(ring==2) hpointNringTwo[0]->Fill(nCrossed[layer]); else if(ring==3) hpointNringThree[0]->Fill(nCrossed[layer]);
	   }
	   ///cout << " No hit in layer " << layer << ": nCrossed = " << nCrossed[layer] << ", nMatched = " << nMatched[layer] << endl;
	 }
       }

       ///cout << " " << endl;

       int nRpc = 0, nRB = 0, nRE = 0;
       int nMatchedRpc = 0, nMatchedRB = 0, nMatchedRE = 0;
       for ( RecToRecHitMap::const_iterator match = refToRecHitMap.begin();
	     match != refToRecHitMap.end(); ++match )
	 {
	   RecHitIter recHitIter = match->first;
	   RecHitIter refHitIter = match->second;

	   const RPCDetId recDetId = static_cast<const RPCDetId>(recHitIter->rpcId());
	   const RPCRoll* recRoll = dynamic_cast<const RPCRoll*>(rpcGeometry->roll(recDetId));

	   const RPCDetId refDetId = static_cast<const RPCDetId>(refHitIter->rpcId());
	   const RPCRoll* refRoll = dynamic_cast<const RPCRoll*>(rpcGeometry->roll(refDetId));

	   RPCGeomServ RPCname(recDetId);
	   const int index     = rpcIndex[RPCname.name()];
	   const int region    = recDetId.region();
	   const int ring      = recDetId.ring();
	   const int sector    = recDetId.sector();
	   const int station   = recDetId.station(); //roll->id().station()
	   const int sublayer  = recDetId.layer();
	   const int subsector = recDetId.subsector();
	   const int roll      = recDetId.roll();
	   const int rawid     = recDetId.rawId();
	   int layer = 0;
	   if (region==0) {
	     layer = station-1 + station*sublayer;
	     if ((station==2 && sublayer==2) || (station==4 && sublayer==1)) layer -= 1;
	   } else if(region==-1) layer = station + 6;
	   else layer = station + 9;

	   int diskTmp = station;
	   if (region==-1) diskTmp += 4; else if (region==1) diskTmp += 7;
	   
           int ringTmp = 0;
           if(region==-1) {
             if(ring==2) {
               if(station==1) ringTmp = 1;       //RE ring -2
               else if(station==2) ringTmp = 2;  //RE ring -2
               else if(station==3) ringTmp = 3;  //RE ring -2
             } else if(ring==3) {
               if(station==1) ringTmp = 7;       //RE ring -3
               else if(station==2) ringTmp = 8;  //RE ring -3
               else if(station==3) ringTmp = 9;  //RE ring -3
             }
           } else if(region==1) {
             if(ring==2) {
               if(station==1) ringTmp = 4;       //RE ring +2
               else if(station==2) ringTmp = 5;  //RE ring +2
               else if(station==3) ringTmp = 6;  //RE ring +2
             } else if(ring==3) {
               if(station==1) ringTmp = 10;      //RE ring +3
               else if(station==2) ringTmp = 11; //RE ring +3
               else if(station==3) ringTmp = 12; //RE ring +3
             }
           }

	   LocalPoint recHitLocal    = recHitIter->localPosition();
	   GlobalPoint rpcGlobalPos  = recRoll->surface().toGlobal(recHitLocal);
	   //const double globalZ      = rpcGlobalPos.z(), globalPerp = rpcGlobalPos.perp();
	   const double globalEta    = rpcGlobalPos.eta(), globalPhi = rpcGlobalPos.phi(); 
	   const double localX       = recHitLocal.x();

	   LocalPoint hitLocal       = refHitIter->localPosition();
	   GlobalPoint hitGlobalPos  = refRoll->surface().toGlobal(hitLocal);
	   //const double hitGlobalZ   = hitGlobalPos.z(), hitGlobalPerp = hitGlobalPos.perp(), hitGlobalEta = hitGlobalPos.eta();
	   const double hitGlobalPhi = hitGlobalPos.phi(); 
	   const double hitLocalX    = hitLocal.x();

	   const double res          = localX - hitLocalX;
	   const double res_phi      = globalPhi - hitGlobalPhi;

           const double localXerr    = recHitIter->localPositionError().xx();
           const double hitLocalXerr = refHitIter->localPositionError().xx();
           const double pull  = (localXerr+hitLocalXerr == 0) ? -999 : res/sqrt(hitLocalXerr); //e.g (Xreco-genX)/sigma(reco)
           //const double pulltmp = (localXerr+hitLocalXerr == 0) ? -999 : res/sqrt(localXerr+hitLocalXerr);
	   ///cout << " Pull in layer " << layer << ": " << pull << " with localX error = " << localXerr << endl ;

           nRpc++; if(region==0) nRB++; else nRE++;

           TString WHEEL = region==0 ? "wheel" : "ring ";
           //cout << " layer " << layer << ": " << WHEEL << " (" << ring << " vs. " << minDxRing[layer] << ") and Residual (" << res << " vs. " << minDx[layer] << "), nMatched = " << nMatched[layer] << endl;
           //if(fabs(res - minDx[layer])>0.5) cout << "Found difference in layer " << layer << endl;
           ///cout << " Hit in layer " << layer << ": " << WHEEL << " (" << ring << " vs. " << minDxRing[layer] << ") and Residual (" << res_phi << " vs. " << minDphi[layer] << "), nCrossed = " << nCrossed[layer] << ", nMatched = " << nMatched[layer] << endl;
           ///if(fabs(res_phi - minDphi[layer])>0.0001) cout << "Found difference in layer " << layer << endl;

	   hresidual[0]->Fill(res); hresidual[diskTmp]->Fill(res); hres[0]->Fill(res); hres[diskTmp]->Fill(res);
	   hpull[0]->Fill(pull); hpull[diskTmp]->Fill(pull);
           if(region==0) {
             hresidualRB[0]->Fill(res); hresidualRB[station]->Fill(res); hresidualRBlayer[layer]->Fill(res);
             hresRB[0]->Fill(res); hresRB[station]->Fill(res); hresRBlayer[layer]->Fill(res);
             hpullRB[0]->Fill(pull); hpullRB[station]->Fill(pull); hpullRBlayer[layer]->Fill(pull);
             htestRB[0]->Fill(res_phi,minDphi[layer]); htestRB[station]->Fill(res_phi,minDphi[layer]); htestRBlayer[layer]->Fill(res_phi,minDphi[layer]);
           } else {
             hresidualRE[0]->Fill(res); hresidualRE[station]->Fill(res); hresidualRing[ring]->Fill(res);
             if(region==-1) hresidualRing[ring-2]->Fill(res); else hresidualRing[ring+2]->Fill(res);

             hresRE[0]->Fill(res); hresRE[station]->Fill(res); hresRing[ring]->Fill(res);
             if(region==-1) hresRing[ring-2]->Fill(res); else hresRing[ring+2]->Fill(res);

             hpullRE[0]->Fill(pull); hpullRE[station]->Fill(pull); hpullRing[ring]->Fill(pull);
             if(region==-1) hpullRing[ring-2]->Fill(pull); else hpullRing[ring+2]->Fill(pull);

             htestRE[0]->Fill(res_phi,minDphi[layer]); htestRE[station]->Fill(res_phi,minDphi[layer]); htestRing[ring]->Fill(res_phi,minDphi[layer]);
             if(region==-1) htestRing[ring-2]->Fill(res_phi,minDphi[layer]); else htestRing[ring+2]->Fill(res_phi,minDphi[layer]);
           }

	   hresRoll[0][index]->Fill(res); hresRoll[1][index]->Fill(res);

	   ///if(fabs(res)>50) cout << " Unmatched hit with large residual = " << res << " in layer " << layer << endl ;
	   //bool associated = region==0 ? fabs(res) < 25 : fabs(res) < 20; //i.e., passed residual cut [cm]
           bool associated = fabs(res) < maxRes;

           if( ! associated ) continue;

           nMatchedRpc++; if(region==0) nMatchedRB++; else nMatchedRE++;

	   heffRPCRecHit[0]->Fill(1); if(region==0) heffRPCRecHit[1]->Fill(1); else heffRPCRecHit[2]->Fill(1);
	   heffEta[1]->Fill(globalEta); heffPhi[1]->Fill(globalPhi); heffEtaPhi[1]->Fill(globalEta,globalPhi);
           heffEta[0]->Fill(globalEta); heffPhi[0]->Fill(globalPhi); heffEtaPhi[0]->Fill(globalEta,globalPhi);

           rpcFound[0][RPCname.name()]++; if(region==0) rpcFound[1][RPCname.name()]++; else rpcFound[2][RPCname.name()]++;
           expected[0][RPCname.name()]++; if(region==0) expected[1][RPCname.name()]++; else expected[2][RPCname.name()]++;

	   for(int i=0; i<2; ++i) {
	     if (region==0) {
	       hlocalXRB[i][0]->Fill(localX); hlocalXRB[i][station]->Fill(localX); hlocalXLayer[i][0]->Fill(localX); hlocalXLayer[i][layer]->Fill(localX);
	     } else {
	       int diskTmp = layer - 6;
	       hlocalXRE[i][0]->Fill(localX); hlocalXRE[i][station]->Fill(localX); htmpXRE[i][0]->Fill(localX); htmpXRE[i][diskTmp]->Fill(localX);
               hlocalXring[i][0]->Fill(localX); hlocalXring[i][ringTmp]->Fill(localX);
               if(ring==2) hlocalXringTwo[i]->Fill(localX); else if(ring==3) hlocalXringThree[i]->Fill(localX);
             }
           }

	 }

       hrechit[2]->Fill(nRpc); hrechitRB[2]->Fill(nRB); hrechitRE[2]->Fill(nRE);
       if(abs(etaMu)<0.8) hrechitRB[3]->Fill(nRpc); else if(abs(etaMu)>1.2) hrechitRE[3]->Fill(nRpc); else hrechitRBRE[3]->Fill(nRpc);

       hrechit[4]->Fill(nMatchedRpc); hrechitRB[4]->Fill(nMatchedRB); hrechitRE[4]->Fill(nMatchedRE);
       if(abs(etaMu)<0.8) hrechitRB[5]->Fill(nMatchedRpc); else if(abs(etaMu)>1.2) hrechitRE[5]->Fill(nMatchedRpc); else hrechitRBRE[5]->Fill(nMatchedRpc);

       //RPC+tracker should have at least one hit in RPC matched to trajectory (no residual cut)
       if(nRpc>0) {
         hmupt[1][0]->Fill(ptMu); hmueta[1][0]->Fill(etaMu); hmuphi[1][0]->Fill(phiMu);
         if(abs(etaMu)<0.8) { hmupt[1][1]->Fill(ptMu); hmuphi[1][1]->Fill(phiMu); }
         else if(abs(etaMu)>1.2) { hmupt[1][3]->Fill(ptMu); hmuphi[1][3]->Fill(phiMu); }
         else { hmupt[1][2]->Fill(ptMu); hmuphi[1][2]->Fill(phiMu); }
       }
	
       //--updated since Nov 22 to get exactly 100% efficiency with zero RPC hit associated  
       for(int k=0; k<9; ++k) {
	 if(k>nMatchedRpc) continue;
	 hpt[1][0][k]->Fill(ptMu); heta[1][0][k]->Fill(etaMu); hphi[1][0][k]->Fill(phiMu);
	 if(abs(etaMu)<0.8) { hpt[1][1][k]->Fill(ptMu); hphi[1][1][k]->Fill(phiMu); }
	 else if(abs(etaMu)>1.2) { hpt[1][3][k]->Fill(ptMu); hphi[1][3][k]->Fill(phiMu); }
	 else { hpt[1][2][k]->Fill(ptMu); hphi[1][2][k]->Fill(phiMu); }
       }

     } //end loop on MuonCollection

     LogVerbatim("TrackAssociator") << "Finally, number of muons selected in the event: " << nMuons << endl ;
     hmuN->Fill(nMuons); //number of muons selected

     for ( int ii=0, nn=rpcAllHitVector.size(); ii<nn; ++ii )
     {
       delete rpcAllHitVector[ii];
     }
     rpcAllHitVector.clear();

   } //EventSelection
}

void TestTrackAssociator::beginJob(){

  theStart=true;
  //theStart=false;

  //bxTree->Branch("hitbx",    &hitbx,    "hitbx/I");
  //bxTree->Branch("hitdx",    &hitdx,    "hitdx/D");
  //bxTree->Branch("detector", &detector, "detector/I");
  //bxTree->Branch("name",     name,      "name/C");

}

void TestTrackAssociator::endJob(){

  honeDivide(heffEta[2],heffEta[1],heffEta[0]);
  honeDivide(heffPhi[2],heffPhi[1],heffPhi[0]);
  htwoDivide(heffEtaPhi[2],heffEtaPhi[1],heffEtaPhi[0]);

  for(int j=0; j<13; ++j) {
    if(j<7) honeDivide(hlocalXLayer[2][j],hlocalXLayer[1][j],hlocalXLayer[0][j]);
    if(j<5) honeDivide(hlocalXRB[2][j],hlocalXRB[1][j],hlocalXRB[0][j]);
    if(j<4) honeDivide(hlocalXRE[2][j],hlocalXRE[1][j],hlocalXRE[0][j]);
    if(j<7) honeDivide(htmpXRE[2][j],htmpXRE[1][j],htmpXRE[0][j]);

    honeDivide(hlocalXring[2][j],hlocalXring[1][j],hlocalXring[0][j]);
    if(j==0) {
      honeDivide(hlocalXringTwo[2],hlocalXringTwo[1],hlocalXringTwo[0]);
      honeDivide(hlocalXringThree[2],hlocalXringThree[1],hlocalXringThree[0]);
    }
  }

  for(int j=0; j<4; ++j) {
    honeDivide(hmupt[2][j],hmupt[1][j],hmupt[0][j]);
    honeDivide(hmueta[2][j],hmueta[1][j],hmueta[0][j]);
    honeDivide(hmuphi[2][j],hmuphi[1][j],hmuphi[0][j]);
    for(int k=0; k<7; ++k) {
      //honeDivide(hpt[2][j][k],hpt[1][j][k],hpt[0][j][k]);
      //honeDivide(heta[2][j][k],heta[1][j][k],heta[0][j][k]);
      //honeDivide(hphi[2][j][k],hphi[1][j][k],hphi[0][j][k]);
      honeDivide(hpt[2][j][k],hpt[1][j][k],hmupt[0][j]);
      honeDivide(heta[2][j][k],heta[1][j][k],hmueta[0][j]);
      honeDivide(hphi[2][j][k],hphi[1][j][k],hmuphi[0][j]);
    }
  }

  // create output file for efficiency
  ofstream f;
  f.open(outputAscii_.c_str());
  //  f.open("effdummy"); 
  if (!f) cout << "Error in Efficiency file open " << endl;

  TFile *file = new TFile("residual.root","RECREATE");
  TTree *rpcTree = new TTree("rpcTree","rpcTree");
  rpcTree->Branch("name",        name,         "name/C");
  rpcTree->Branch("detector",    &detector,    "detector/I");
  rpcTree->Branch("denominator", &denominator, "denominator/I");
  rpcTree->Branch("numerator",   &numerator,   "numerator/I");
  rpcTree->Branch("efficiency",  &efficiency,  "efficiency/D");
  rpcTree->Branch("error",       &error,       "error/D");

  map<string, int>::const_iterator iter;
  //int denominator, numerator;
  //float efficiency, error;
  denominator = numerator = efficiency = error = detector = -999;
  for (int i=0; i<3; ++i) {
    for (iter=expected[i].begin(); iter != expected[i].end(); ++iter) {
      string nameString = iter->first;
      sscanf(nameString.c_str(),"%s",name);

      TString nameChamber = nameString;
      if(nameChamber.Contains("RB")) detector = 1;
      else if(nameChamber.Contains("RE")) detector = 2;

      denominator = iter->second;
      numerator = rpcFound[i][iter->first];
      efficiency=0.; error = 0.;
      if (denominator!=0) {
	efficiency = float(numerator)/float(denominator);
	error = sqrt(efficiency*(1-efficiency)/(float)denominator);
	heffRPC[i]->Fill(efficiency); herrRPC[i]->Fill(error); heffErr[i]->Fill(efficiency,error);

        if(i==0) {
	  f << name << " " << denominator << " " << numerator << " "  << efficiency << " " << error << endl;
	  rpcTree->Fill();
	}

      }
    }
  }
  f.close();

  //double mean, rms;
  //int nameId;
  //string name;

  TNtuple *ntuple = new TNtuple("ntuple","residuals","mean:rms:entries:nameId");
  TTree *T = new TTree("T","tree");
  T->Branch("mean",    &mean,    "mean/D");
  T->Branch("rms",     &rms,     "rms/D");
  T->Branch("mean2",   &mean2,   "mean2/D");
  T->Branch("rms2",    &rms2,    "rms2/D");
  T->Branch("entries", &entries, "entries/I");
  T->Branch("nameId",  &nameId,  "nameId/I");
  T->Branch("name",    name,     "name/C");

  map<string, int>::const_iterator iterNew;
  for (int i=0; i<3; ++i) {
    //for (int j=0; j<2316; ++j) { nameId = j;
    for (iterNew=rpcIndex.begin(); iterNew != rpcIndex.end(); ++iterNew) {
      string nameString = iterNew->first; //name = nameString.c_str();
      sscanf(nameString.c_str(),"%s",name);
      nameId = iterNew->second;

      if(i<2) { mean = hresRoll[i][nameId]->GetMean(); rms = hresRoll[i][nameId]->GetRMS(); entries = hresRoll[i][nameId]->GetEntries(); }
      else if(i==1) { mean2 = hresRoll[i][nameId]->GetMean(); rms2 = hresRoll[i][nameId]->GetRMS(); }
      else { mean = hresRoll[0][nameId]->GetMean(); rms = hresRoll[0][nameId]->GetRMS(); entries = hresRoll[0][nameId]->GetEntries(); }
      hresMean[i]->Fill(mean); hresRms[i]->Fill(rms); hresEntries[i]->Fill(entries);
      hresMeanRms[i]->Fill(rms,mean); hresMeanIndex[i]->Fill(nameId,mean); hresRmsIndex[i]->Fill(nameId,rms);
      if (i==0) {
        //cout << nameId << " " << name << " " << mean << " " << rms << endl;
        T->Fill(); ntuple->Fill(mean,rms,entries,nameId);
      }
    }
  }
  file->Write();
  file->Close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(TestTrackAssociator);
