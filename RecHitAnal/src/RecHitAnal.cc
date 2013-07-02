// -*- C++ -*-                                                                                                                   
//                                                                                                                               
// Package:    RecHitAnal                                                                                                            
// Class:      RecHitAnal                                                                                                            
//                                                                                                                               
/**\class RecHitAnal RecHitAnal.cc RPCAnalysis/RecHitAnal/src/RecHitAnal.cc                                                                      
                                                                                                                                 
 Description: [one line class summary]                                                                                           
                                                                                                                                 
 Implementation:                                                                                                                 
     [Notes on implementation]                                                                                                   
*/
//                                                                                                                               
// Original Author:  Hyunkwan Seo,588 R-009,+41227678393,                                                                        
//         Created:  Fri Jun 17 17:45:39 CEST 2011<<<<<<< RecHitAnal.cc
// $Id: RecHitAnal.cc,v 1.11 2012/12/14 16:43:20 hkseo Exp $=======
// $Id: RecHitAnal.cc,v 1.11 2012/12/14 16:43:20 hkseo Exp $>>>>>>> 1.4
//                                                                                                                               
//    

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/SegmentsTrackAssociator.h"
#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"

#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCChamber.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DQM/RPCMonitorDigi/interface/utils.h"
#include "DQMOffline/Muon/src/SegmentTrackAnalyzer.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// global trigger stuff
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerEvmReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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


using namespace edm; 
using namespace std;
using namespace reco;

static const int kMax = 10;
static const int kMaxSegm = 20;
static const int kMaxMu = 10000;

//
// class declaration
//
 
class RecHitAnal : public edm::EDAnalyzer {
public:
  explicit RecHitAnal(const edm::ParameterSet&);
  ~RecHitAnal();
 
 
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  // ----------member data ---------------------------
  int fBX;
  int fOrbit;

  //Tree for each event
  TTree *t0; 
  Int_t Run;
  Int_t eventNumber;
  Int_t totalEvents;
  Float_t angle;
  Int_t nTracks;
  Int_t nMuons;
  Int_t nCSCseg;
  Int_t nDTseg;
  Float_t dxy[kMaxMu]; 
  Float_t dz_[kMaxMu];
  Float_t mupt[kMaxMu];
  Float_t mup[kMaxMu];
  Float_t mupx[kMaxMu];
  Float_t mupy[kMaxMu];
  Float_t mupz[kMaxMu];
  Float_t muphi[kMaxMu];
  Float_t mueta[kMaxMu];
  Float_t muet[kMaxMu];
  Float_t mue[kMaxMu];
  Int_t muchg[kMaxMu];

  bool isTrackerMu[kMaxMu];
  bool isStaMu[kMaxMu];
  bool isGlobalMu[kMaxMu];
  bool isPFMu[kMaxMu];
  bool isMatchToNoRPC[kMaxMu];

  // Tree for each muon
  TTree *t1; 
  Int_t hitsFromRpc, hitsFromDt, hitsFromCsc;
  Int_t hitsFromRpcSTA, hitsFromDtSTA, hitsFromCscSTA;
  Int_t nValidMuonRPCHits, rpcStationsWithValidHits;
  Int_t nValidHits, nValidMuonHits, nValidTrackerHits, nValidPixelHits, nTrackerLayers;
  int segmFromDt, segmFromCsc;
  Int_t nMatch; //new

  Int_t nValidMuonHitsNoRPC, nValidTrackerHitsNoRPC; //new
  Int_t segmFromDtNoRPC, segmFromCscNoRPC;//new
  Int_t nMatchNoRPC, nMatchedSegmNoRPC;//new

  Float_t eta, phi, pt;
  Float_t eta_STA, phi_STA, pt_STA;
  Float_t normChi2;
  Double_t d0,dz;
  Bool_t matchToNoRPC;
  Bool_t isSTA, isPF;

  //for GLB
  Int_t region[kMax];
  Int_t ring[kMax];
  Int_t sector[kMax];
  Int_t subsector[kMax];
  Int_t station[kMax];  // 1-4 for barrel, disk 1-3 for endcaps 
  Int_t sublayer[kMax]; // 1 (in) or 2 (out) just for barrel RB1 and RB2  
  Float_t recX[kMax];
  Int_t hitsFromSegmDt[kMaxSegm];
  Int_t hitsFromSegmCsc[kMaxSegm];
  Float_t chi2FromSegmDt[kMaxSegm];
  Float_t chi2FromSegmCsc[kMaxSegm];
  Float_t ndofFromSegmDt[kMaxSegm];
  Float_t ndofFromSegmCsc[kMaxSegm];

  Int_t nMatchedSegm;


  //for STA
  Int_t regionSTA[kMax];
  Int_t ringSTA[kMax];
  Int_t sectorSTA[kMax];
  Int_t subsectorSTA[kMax];
  Int_t stationSTA[kMax];  // 1-4 for barrel, disk 1-3 for endcaps 
  Int_t sublayerSTA[kMax]; // 1 (in) or 2 (out) just for barrel RB1 and RB2  
  Float_t recXSTA[kMax];


  //for TRK Muons
  TTree *t3; 
  Int_t hitsFromRpc_track, hitsFromDt_track, hitsFromCsc_track;
  Int_t segmFromDt_track, segmFromCsc_track;
  Float_t pt_track, eta_track, phi_track;
  
  Int_t region_track[kMax];
  Int_t ring_track[kMax];
  Int_t sector_track[kMax];
  Int_t subsector_track[kMax];
  Int_t station_track[kMax];  // 1-4 for barrel, disk 1-3 for endcaps 
  Int_t sublayer_track[kMax]; // 1 (in) or 2 (out) just for barrel RB1 and RB2  
  Float_t recX_track[kMax];
  Int_t hitsFromSegmDt_track[kMaxSegm];
  Int_t hitsFromSegmCsc_track[kMaxSegm];
  Float_t chi2FromSegmDt_track[kMaxSegm];
  Float_t chi2FromSegmCsc_track[kMaxSegm];
  Float_t ndofFromSegmDt_track[kMaxSegm];
  Float_t ndofFromSegmCsc_track[kMaxSegm];

  bool Debug_;

  InputTag RPCRecHits_;
  InputTag DT4DSegments_;
  InputTag CSCSegments_;
  InputTag muons_;
  InputTag globalMuonsNoRPC_;
  InputTag StandAloneMuons_;
  InputTag genTracksLabel_;

  SegmentsTrackAssociator* theSegmentsAssociator;


  // Tree for each CSC segment
  TTree *t2; 
  Int_t endcap, stationCSC, ringCSC, chamberCSC;
  Float_t globalX, globalY;

  //Service<TFileService> fs;
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecHitAnal::RecHitAnal(const edm::ParameterSet& cfg)

{
  
   //now do what ever initialization is needed
  Service<TFileService> fs;

  const ParameterSet SegmentsTrackAssociatorParameters = cfg.getParameter<ParameterSet>("SegmentsTrackAssociatorParameters");
  theSegmentsAssociator = new SegmentsTrackAssociator(SegmentsTrackAssociatorParameters);
  //theSegmentsAssociator = new SegmentsTrackAssociator("SegmentsTrackAssociatorParameters");
  
  //debug
  Debug_ = cfg.getUntrackedParameter<bool>("Debug",false);
  
  //RPC hits and digis 
  RPCRecHits_ = cfg.getParameter<InputTag>("RPCRecHits"); 
  DT4DSegments_ = cfg.getParameter<InputTag>("DT4DSegments");
  CSCSegments_ = cfg.getParameter<InputTag>("CSCSegments");
  muons_ = cfg.getParameter<InputTag>("muons");
  globalMuonsNoRPC_ = cfg.getParameter<InputTag>("globalMuonsNoRPC");
  StandAloneMuons_ = cfg.getParameter<InputTag>("StandAloneMuons");
  genTracksLabel_ = cfg.getParameter<InputTag>("generalTracks");


  //create the Tree and branches
  t0 = fs->make<TTree>("Events","");
  t0->Branch("Run",         &Run,              "Run/I");
  t0->Branch("eventNumber", &eventNumber,      "eventNumber/I");
  t0->Branch("totalEvents", &totalEvents,      "totalEvents/I");
  t0->Branch("angle",       &angle,            "angle/F");
  t0->Branch("nTracks",     &nTracks,          "nTracks/I");
  t0->Branch("nCSCseg",     &nCSCseg,          "nCSCseg/I");
  t0->Branch("nDTseg",      &nDTseg,           "nDTseg/I");
  t0->Branch("NrecoMuon",   &nMuons,           "NrecoMuon/I");
  t0->Branch("dxy",         dxy,               "dxy[NrecoMuon]/F");
  t0->Branch("dz",          dz_,               "dz[NrecoMuon]/F");
  t0->Branch("recoMuonPt",  mupt,              "recoMuonPt[NrecoMuon]/F");
  t0->Branch("recoMuonP",   mup,               "recoMuonP[NrecoMuon]/F");
  t0->Branch("recoMuonPx",  mupx,              "recoMuonPx[NrecoMuon]/F");
  t0->Branch("recoMuonPy",  mupy,              "recoMuonPy[NrecoMuon]/F");
  t0->Branch("recoMuonPz",  mupz,              "recoMuonPz[NrecoMuon]/F");
  t0->Branch("recoMuonPhi", muphi,             "recoMuonPhi[NrecoMuon]/F");
  t0->Branch("recoMuonEta", mueta,             "recoMuonEta[NrecoMuon]/F");
  t0->Branch("recoMuonEt",  muet,              "recoMuonEt[NrecoMuon]/F");
  t0->Branch("recoMuonE",   mue,               "recoMuonE[NrecoMuon]/F");
  t0->Branch("isTrackerMu", isTrackerMu,       "isTrackerMu[NrecoMuon]/O");
  t0->Branch("isStaMu",     isStaMu,           "isStaMu[NrecoMuon]/O");
  t0->Branch("isGlobalMu",  isGlobalMu,        "isGlobalMu[NrecoMuon]/O");
  t0->Branch("isPFMu",      isPFMu,            "isPFMu[NrecoMuon]/O");
  t0->Branch("isMatchToNoRPC",  isMatchToNoRPC,  "isMatchToNoRPC[NrecoMuon]/O");


  t1 = fs->make<TTree>("muon","");
  t1->Branch("Run",         &Run,              "Run/I");
  t1->Branch("eventNumber", &eventNumber,      "eventNumber/I");
  t1->Branch("angle",       &angle,            "angle/F");
  t1->Branch("dxy",         &d0,               "dxy/D");
  t1->Branch("dz",          &dz,               "dz/D");
  t1->Branch("matchToNoRPC",&matchToNoRPC,     "matchToNoRPC/O");
  t1->Branch("isStaMu",     &isSTA,            "isStaMu/O");
  t1->Branch("isPFMu",      &isPF,             "isPFMu/O");
  t1->Branch("nTracks",     &nTracks,          "nTracks/I");
  t1->Branch("nMuons",      &nMuons,           "nMuons/I");
  t1->Branch("nRpcHit",     &hitsFromRpc,      "nRpcHit/I");


  t1->Branch("nValidHits",              &nValidHits,               "nValidHits/I");
  t1->Branch("nValidMuonHits",          &nValidMuonHits,           "nValidMuonHits/I");
  t1->Branch("nValidTrackerHits",       &nValidTrackerHits,        "nValidTrackerHits/I");
  t1->Branch("nValidPixelHits",         &nValidPixelHits,          "nValidPixelHits/I");
  t1->Branch("nTrackerLayers",          &nTrackerLayers,           "nTrackerLayers/I");
  t1->Branch("nValidMuonRpcHit",        &nValidMuonRPCHits,        "nValidMuonRpcHit/I");
  t1->Branch("nMatch",            &nMatch,           "nMatch/I");
  t1->Branch("nMatchedSegm",      &nMatchedSegm,     "nMatchedSegm/I");

  t1->Branch("nValidMuonHitsNoRPC",     &nValidMuonHitsNoRPC,      "nValidMuonHitsNoRPC/I"); //new
  t1->Branch("nValidTrackerHitsNoRPC",  &nValidTrackerHitsNoRPC,   "nValidTrackerHitsNoRPC/I"); //new

  t1->Branch("nMatchNoRPC",       &nMatchNoRPC,           "nMatchNoRPC/I"); //new
  t1->Branch("nMatchedSegmNoRPC", &nMatchedSegmNoRPC,     "nMatchedSegmNoRPC/I"); //new
  t1->Branch("nDtSegmNoRPC",      &segmFromDtNoRPC,       "nDtSegmNoRPC/I"); //new
  t1->Branch("nCscSegmNoRPC",     &segmFromCscNoRPC,      "nCscSegmNoRPC/I"); //new


  t1->Branch("rpcStationsWithHits",     &rpcStationsWithValidHits, "rpcStationsWithHits/I");
  t1->Branch("nDtHit",      &hitsFromDt,       "nDtHit/I");
  t1->Branch("nCscHit",     &hitsFromCsc,      "nCscHit/I");
  t1->Branch("nDtSegm",     &segmFromDt,       "nDtSegm/I");
  t1->Branch("nCscSegm",    &segmFromCsc,      "nCscSegm/I");
  t1->Branch("eta",         &eta,              "eta/F");
  t1->Branch("phi",         &phi,              "phi/F");
  t1->Branch("pt",          &pt,               "pt/F");
  t1->Branch("normChi2",    &normChi2,         "normChi2/F");
  t1->Branch("nRpcHitSTA",  &hitsFromRpcSTA,   "nRpcHitSTA/I");
  t1->Branch("nDtHitSTA",   &hitsFromDtSTA,    "nDtHitSTA/I");
  t1->Branch("nCscHitSTA",  &hitsFromCscSTA,   "nCscHitSTA/I");
  t1->Branch("etaSTA",      &eta_STA,          "eta/F");
  t1->Branch("phiSTA",      &phi_STA,          "phi/F");
  t1->Branch("ptSTA",       &pt_STA,           "pt/F");

  // For global muons RPCHits
  t1->Branch("region",          region,            "region[nRpcHit]/I");
  t1->Branch("ring",            ring,              "ring[nRpcHit]/I");
  t1->Branch("sector",          sector,            "sector[nRpcHit]/I");
  t1->Branch("subsector",       subsector,         "subsector[nRpcHit]/I");
  t1->Branch("station",         station,           "station[nRpcHit]/I");
  t1->Branch("sublayer",        sublayer,          "sublayer[nRpcHit]/I");
  t1->Branch("localX",          recX,              "localX[nRpcHit]/F");
  t1->Branch("hitsFromSegmDt",  hitsFromSegmDt,    "hitsFromSegmDt[nDtSegm]/I");
  t1->Branch("hitsFromSegmCsc", hitsFromSegmCsc,   "hitsFromSegmCsc[nCscSegm]/I");
  t1->Branch("chi2FromSegmDt",  chi2FromSegmDt,    "chi2FromSegmDt[nDtSegm]/F");
  t1->Branch("chi2FromSegmCsc", chi2FromSegmCsc,   "chi2FromSegmCsc[nCscSegm]/F");
  t1->Branch("ndofFromSegmDt",  ndofFromSegmDt,    "ndofFromSegmDt[nDtSegm]/F");
  t1->Branch("ndofFromSegmCsc", ndofFromSegmCsc,   "ndofFromSegmCsc[nCscSegm]/F");

  // For standalone muons RPCHits
  t1->Branch("regionSTA",   regionSTA,         "regionSTA[nRpcHitSTA]/I");
  t1->Branch("ringSTA",     ringSTA,           "ringSTA[nRpcHitSTA]/I");
  t1->Branch("sectorSTA",   sectorSTA,         "sectorSTA[nRpcHitSTA]/I");
  t1->Branch("subsectorSTA",subsectorSTA,      "subsectorSTA[nRpcHitSTA]/I");
  t1->Branch("stationSTA",  stationSTA,        "stationSTA[nRpcHitSTA]/I");
  t1->Branch("sublayerSTA", sublayerSTA,       "sublayerSTA[nRpcHitSTA]/I");
  t1->Branch("localXSTA",   recXSTA,           "localXSTA[nRpcHitSTA]/F");

  t2 = fs->make<TTree>("CSCSegm","");
  t2->Branch("Run",         &Run,              "Run/I");
  t2->Branch("eventNumber", &eventNumber,      "eventNumber/I");
  t2->Branch("totalEvents", &totalEvents,      "totalEvents/I");
  t2->Branch("endcap",      &endcap,           "endcap/I");
  t2->Branch("station",     &stationCSC,       "station/I");
  t2->Branch("ring",        &ringCSC,          "ring/I");
  t2->Branch("chamber",     &chamberCSC,       "chamber/I");
  t2->Branch("globalX",     &globalX,          "globalX/F");
  t2->Branch("globalY",     &globalY,          "globalY/F");

  t3 = fs->make<TTree>("TrackerMuon","");
  t3->Branch("nRpcHit_track",  &hitsFromRpc_track,   "nRpcHit_track/I");
  t3->Branch("nDtHit_track",   &hitsFromDt_track,    "nDtHit_track/I");
  t3->Branch("nCscHit_track",  &hitsFromCsc_track,   "nCscHit_track/I");
  t3->Branch("eta_track",      &eta_track,        "eta_track/F");
  t3->Branch("phi_track",      &phi_track,        "phi_track/F");
  t3->Branch("pt_track",       &pt_track,         "pt_track/F");
  t3->Branch("region_track",   region_track,         "region_track[nRpcHit_track]/I");
  t3->Branch("ring_track",     ring_track,           "ring_track[nRpcHit_track]/I");
  t3->Branch("sector_track",   sector_track,         "sector_track[nRpcHit_track]/I");
  t3->Branch("subsector_track",subsector_track,      "subsector_track[nRpcHit_track]/I");
  t3->Branch("station_track",  station_track,        "station_track[nRpcHit_track]/I");
  t3->Branch("sublayer_track", sublayer_track,       "sublayer_track[nRpcHit_track]/I");
  t3->Branch("localX_track",   recX_track,           "localX_track[nRpcHit_track]/F");

  t3->Branch("nDtSegm_track",     &segmFromDt_track,       "nDtSegm_track/I");
  t3->Branch("nCscSegm_track",    &segmFromCsc_track,      "nCscSegm_track/I");

  t3->Branch("hitsFromSegmDt_track",  hitsFromSegmDt_track,    "hitsFromSegmDt_track[nDtSegm_track]/I");
  t3->Branch("hitsFromSegmCsc_track", hitsFromSegmCsc_track,   "hitsFromSegmCsc_track[nCscSegm_track]/I");
  t3->Branch("chi2FromSegmDt_track",  chi2FromSegmDt_track,    "chi2FromSegmDt_track[nDtSegm_track]/F");
  t3->Branch("chi2FromSegmCsc_track", chi2FromSegmCsc_track,   "chi2FromSegmCsc_track[nCscSegm_track]/F");
  t3->Branch("ndofFromSegmDt_track",  ndofFromSegmDt_track,    "ndofFromSegmDt_track[nDtSegm_track]/F");
  t3->Branch("ndofFromSegmCsc_track", ndofFromSegmCsc_track,   "ndofFromSegmCsc_track[nCscSegm_track]/F");

}


RecHitAnal::~RecHitAnal()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
 
}

//
// member functions
//


// ------------ method called to for each event  ------------
void RecHitAnal::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  
  using namespace edm;
  if(Debug_) cout<< " mark 0 "<<endl;

  // initialization for every event


  // the event
  Run = event.id().run();
  eventNumber = event.id().event();
  totalEvents = 1;
  fBX        = event.bunchCrossing();  // LHC bunch crossing
  fOrbit     = event.orbitNumber();
  //cout << "BX " << fBX << "   orbit " << fOrbit << endl;
  //cout <<Run<<" "<<eventNumber<<endl;

 
  // Get the DT 4D segment collection from the event
  Handle<DTRecSegment4DCollection> allDT4DSegments;
  event.getByLabel(DT4DSegments_, allDT4DSegments);

  // Get the csc segment collection from the event
  Handle<CSCSegmentCollection> allCSCSegments;
  event.getByLabel(CSCSegments_, allCSCSegments);
 
  //RPC hits
  Handle<RPCRecHitCollection > rpcRecHits;
  event.getByLabel(RPCRecHits_,rpcRecHits);

  //Access to the beam spot
  reco::BeamSpot beamSpot;
  Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel("offlineBeamSpot", beamSpotHandle);
  
  if ( beamSpotHandle.isValid() )
    {
      beamSpot = *beamSpotHandle;
      
    } else
    {
      cout << "No beam spot available from EventSetup \n" << endl;
    }
  
  //Input standalone muons
  Handle<View<Track> > StandAloneMuons;
  event.getByLabel(StandAloneMuons_,StandAloneMuons);

  //Input muons
  Handle< View<Muon> > muons;      
  event.getByLabel(muons_,muons);

  Handle< View<Muon> > globalMuonsNoRPC;
  event.getByLabel(globalMuonsNoRPC_,globalMuonsNoRPC);


  //Input general tracks
  Handle<reco::TrackCollection> generalTracks;
  event.getByLabel(genTracksLabel_, generalTracks);


  //Geometry, magnetic field and tracks
  ESHandle<MagneticField> MagneticField;
  eventSetup.get<IdealMagneticFieldRecord>().get(MagneticField);
  
  ESHandle<DTGeometry> dtGeometry;
  eventSetup.get<MuonGeometryRecord>().get(dtGeometry);

  ESHandle<CSCGeometry> cscGeometry;
  eventSetup.get<MuonGeometryRecord>().get(cscGeometry);

  ESHandle<RPCGeometry> rpcGeometry;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeometry);

  // Retrieve RecHits from the event
  Handle<RPCRecHitCollection> recHitHandle;
  event.getByLabel("rpcRecHits", recHitHandle);

  typedef RPCRecHitCollection::const_iterator RecHitIter;
  ///  int nStandAloneMuons = StandAloneMuons->size();

  nCSCseg = allCSCSegments->size();
  nDTseg = allDT4DSegments->size();
  //int nRPCRecHits = rpcRecHits->size();
  nMuons = muons->size();
  int nMuonsNoRPC = globalMuonsNoRPC->size();
  //int nStaMuons = StandAloneMuons->size();
  nTracks = generalTracks->size();

  if (nMuons<1) return;   // for now ask for at least one muon


  //*****************************************************************
  //*  Calculate 3D angle for muons in order to remove cosmic muons *
  //*****************************************************************

  vector<TLorentzVector> zmumu;
  angle = -1;

  for (int muidx=0; muidx< nMuons; muidx++) {  // loop on all muons

    TrackRef glbOfGlobalRef = (*muons)[muidx].globalTrack();
    float eta_ = ((*muons)[muidx]).eta();
    //float p_ = ((*muons)[muidx]).p();
    float px_ = ((*muons)[muidx]).px();
    float py_ = ((*muons)[muidx]).py();
    float pz_ = ((*muons)[muidx]).pz();
    float energy_ = ((*muons)[muidx]).energy();
    //int charge_ = ((*muons)[muidx]).charge();
    
    if ( (*muons)[muidx].isGlobalMuon() && abs(eta_)<2.4) {
      math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
      dxy[muidx] = -1.* glbOfGlobalRef->dxy(point);
      zmumu.push_back( TLorentzVector(px_,py_,pz_,energy_) );
    }

  }


  for ( std::vector<TLorentzVector>::iterator mu1 = zmumu.begin(); mu1 != zmumu.end(); ++mu1) {
    for ( std::vector<TLorentzVector>::iterator mu2 = mu1+1; mu2 != zmumu.end(); ++mu2) {
      //double r = mu1->P()*mu2->P();
      //double s = (*mu1).Dot(*mu2); // s = -(*mu1).Dot(*mu2), x = s/r
      //double x = -s/r;
      //double angle3dGLB = acos(x);
      float angleTmp = (*mu1).Angle((*mu2).Vect());
      if (angleTmp > angle) angle = angleTmp;
      //if (angle >= 3.05) return;
    }
  }
  ///////////////////////////////////////////////////////////////////  
  // Finished calculation for 3D angle for muons 
  ///////////////////////////////////////////////////////////////////  

  /*
  cout << "Run " << Run << "  event " << eventNumber << endl;
  cout << " n CSC segments " << nCSCseg << endl;
  cout << " n DT segments " << nDTseg << endl;
  cout << " n RPC hits " << nRPCRecHits << endl;
  cout << " n Muons " << nMuons << endl;
  cout << " n Sta Muons " << nStaMuons << endl;
  */

  vector<GlobalPoint> rpcGlobalPoints;
  for (int muidx=0; muidx< nMuons; muidx++) {  // loop on all muons
    if(Debug_) cout<< " mark 1 "<<endl;
    isTrackerMu[muidx] = false;
    isStaMu[muidx] = false;
    isGlobalMu[muidx] = false;
    isMatchToNoRPC[muidx] = false;
    dxy[muidx] = 0;

    // vector of RpcRecHits associated to sta muon
    //    vector<RPCRecHit> RPCHistOnSta_List;

    mupt[muidx] = ((*muons)[muidx]).pt();
    mup[muidx] = ((*muons)[muidx]).p();
    mupx[muidx] = ((*muons)[muidx]).px();
    mupy[muidx] = ((*muons)[muidx]).py();
    mupz[muidx] = ((*muons)[muidx]).pz();
    muet[muidx] = ((*muons)[muidx]).et();
    mue[muidx] = ((*muons)[muidx]).energy();
    mueta[muidx] = ((*muons)[muidx]).eta();
    muphi[muidx] = ((*muons)[muidx]).phi();
    muchg[muidx] = ((*muons)[muidx]).charge();
    
    eta = ((*muons)[muidx]).eta();
    phi = ((*muons)[muidx]).phi();
    pt = ((*muons)[muidx]).pt();
    //float p = ((*muons)[muidx]).p();
    //float charge = ((*muons)[muidx]).charge();
    nMatch = ((*muons)[muidx]).numberOfMatches();


    //---------- Get Matched Segments --------------
    //     static const int DT= 1;
    //     static const int CSC=2;
    //     static const int RPC=3;
    //nMatchedSegm = 0;

    int csc=0, dt=0;
    for (int station=1; station<5; station++){
      dt += ((*muons)[muidx]).numberOfSegments(station, 1);
      csc += ((*muons)[muidx]).numberOfSegments(station, 2);
    }
    nMatchedSegm = dt + csc;

    //((*muons)[muidx]).matchedSectors(ringMatched, stationMatched, sectorMatched, detIdMatched);

    nMatchNoRPC = -1;
    nMatchedSegmNoRPC = -1;
    nValidMuonHitsNoRPC = -1;
    nValidTrackerHitsNoRPC = -1;
    segmFromDtNoRPC = -1;
    segmFromCscNoRPC = -1;


    TrackRef glbOfGlobalRef = (*muons)[muidx].globalTrack();
    if ( (*muons)[muidx].isTrackerMuon() ) isTrackerMu[muidx] = true;
    if ( (*muons)[muidx].isStandAloneMuon() ) isStaMu[muidx] = true;
    if ( (*muons)[muidx].isGlobalMuon() ) isGlobalMu[muidx] = true;
    if ( (*muons)[muidx].isGlobalMuon() && abs(eta)<2.4) {
      math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
      dxy[muidx] = -1.* glbOfGlobalRef->dxy(point);
      dz_[muidx] = glbOfGlobalRef->dz(point);
    }

    //cout << "   with RPC "<<(*muons)[muidx].isGlobalMuon() << " "<<(*muons)[muidx].isPFMuon() << endl;
    //cout << "without RPC "<<(*globalMuonsNoRPC)[muidx].isGlobalMuon() << " "<<(*globalMuonsNoRPC)[muidx].isPFMuon() << endl;


    //check if muon is matched to mNoRPC
    for(int i=0; i<nMuonsNoRPC; i++) {
      TrackRef glbOfGlobalRefNoRPC = (*globalMuonsNoRPC)[i].globalTrack();
      TLorentzVector t1(mupx[muidx], mupy[muidx], mupz[muidx], mup[muidx]);
      TLorentzVector t2( ((*globalMuonsNoRPC)[i]).px(), ((*globalMuonsNoRPC)[i]).py(), ((*globalMuonsNoRPC)[i]).pz(), ((*globalMuonsNoRPC)[i]).p() );
      if (t2.DeltaR(t1)<0.01) {
        isMatchToNoRPC[muidx] = true;
        nMatchNoRPC = ((*globalMuonsNoRPC)[i]).numberOfMatches();
        nValidMuonHitsNoRPC = glbOfGlobalRefNoRPC->hitPattern().numberOfValidMuonHits();
        nValidTrackerHitsNoRPC = glbOfGlobalRefNoRPC->hitPattern().numberOfValidTrackerHits();
	
        //     static const int DT= 1;
        //     static const int CSC=2;
        //     static const int RPC=3;

        int csc=0, dt=0;
        for (int station=1; station<5; station++){
          dt += ((*globalMuonsNoRPC)[i]).numberOfSegments(station, 1);
          csc += ((*globalMuonsNoRPC)[i]).numberOfSegments(station, 2);
        }
        nMatchedSegmNoRPC = dt + csc;

        // get segments from track
	MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(event, eventSetup, (reco::Track)*glbOfGlobalRefNoRPC);

        // segment counters
        segmFromDtNoRPC=0;
        segmFromCscNoRPC=0;

        for (MuonTransientTrackingRecHit::MuonRecHitContainer::const_iterator segment=segments.begin(); segment!=segments.end(); segment++) {
          DetId id = (*segment)->geographicalId();
          //int hitsFromSegmDt_=0;
          //int hitsFromSegmCsc_=0;

          // hits from DT segments
          if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
            //const DTRecSegment4D *seg4D = dynamic_cast<const DTRecSegment4D*>((*segment)->hit());
            /*
	      if((*seg4D).hasPhi())
              hitsFromSegmDt_+=(*seg4D).phiSegment()->specificRecHits().size();
	      if((*seg4D).hasZed())
              hitsFromSegmDt_+=(*seg4D).zSegment()->specificRecHits().size();
	      hitsFromSegmDt[segmFromDt] = hitsFromSegmDt_;
	      chi2FromSegmDt[segmFromDt] = (*seg4D).chi2();
	      ndofFromSegmDt[segmFromDt] = (*seg4D).degreesOfFreedom();
            */
            ++segmFromDtNoRPC;
          }
          // hits from CSC segments
          if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) {
            /*
	      hitsFromSegmCsc_+=(*segment)->recHits().size();
	      hitsFromSegmCsc[segmFromCsc] = hitsFromSegmCsc_;
	      chi2FromSegmCsc[segmFromCsc] = (*segment)->chi2();
	      ndofFromSegmCsc[segmFromCsc] = (*segment)->degreesOfFreedom();
            */
            segmFromCscNoRPC++;
          }
        }    // end for muonTransientTrackingRecHit
      }  //end if (t2.DeltaR(t1)<0.01)
      
    }
    matchToNoRPC = isMatchToNoRPC[muidx];
    isSTA = isStaMu[muidx];
    
    //if (isGlobalMu[muidx] && isStaMu[muidx] && abs(eta)<1.6) {
    if (isGlobalMu[muidx] && abs(eta)<1.8) {
      if(Debug_) cout<< " mark 2 "<<endl;

      nValidMuonRPCHits = glbOfGlobalRef->hitPattern().numberOfValidMuonRPCHits();
      nValidHits = glbOfGlobalRef->hitPattern().numberOfValidHits();
      nValidMuonHits = glbOfGlobalRef->hitPattern().numberOfValidMuonHits();
      nValidTrackerHits = glbOfGlobalRef->hitPattern().numberOfValidTrackerHits();
      nValidPixelHits = glbOfGlobalRef->hitPattern().numberOfValidPixelHits();
      nTrackerLayers = glbOfGlobalRef->hitPattern().trackerLayersWithMeasurement();
      rpcStationsWithValidHits = glbOfGlobalRef->hitPattern().rpcStationsWithValidHits();
      normChi2 = glbOfGlobalRef->normalizedChi2();

      /*
      float pt_track = trkOfGlobalRef->pt();
      float p_track = trkOfGlobalRef->p();
      float eta_track = trkOfGlobalRef->eta();
      float phi_track = trkOfGlobalRef->phi();
      float chi2_track = trkOfGlobalRef->chi2();
      float normChi2_track = trkOfGlobalRef->normalizedChi2();
      */

      d0 = dxy[muidx];
      dz = dz_[muidx];

      //cout << "dxy: "<<glbOfGlobalRef->dxy()<<" "<<d0<< endl;

      // temporary part
      //      pt = pt_track;
      //      p = p_track;
      //      eta = eta_track;

      // get segments from track
      MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(event, eventSetup, (reco::Track)*glbOfGlobalRef);
      // hit counters
      hitsFromDt=0;
      hitsFromCsc=0;
      hitsFromRpc=0;
      //int hitsFromTk=0;
      int hitsFromTrack=0;
      int validHitsFromTrack=0;

      // segment counters
      segmFromDt=0;
      segmFromCsc=0;

      for (MuonTransientTrackingRecHit::MuonRecHitContainer::const_iterator segment=segments.begin(); segment!=segments.end(); segment++) {
        DetId id = (*segment)->geographicalId();
	int hitsFromSegmDt_=0;
	int hitsFromSegmCsc_=0;	

        // hits from DT segments
        if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
          const DTRecSegment4D *seg4D = dynamic_cast<const DTRecSegment4D*>((*segment)->hit());
          if((*seg4D).hasPhi())
            hitsFromSegmDt_+=(*seg4D).phiSegment()->specificRecHits().size();
          if((*seg4D).hasZed())
            hitsFromSegmDt_+=(*seg4D).zSegment()->specificRecHits().size();
	  hitsFromSegmDt[segmFromDt] = hitsFromSegmDt_;
	  chi2FromSegmDt[segmFromDt] = (*seg4D).chi2(); 
	  ndofFromSegmDt[segmFromDt] = (*seg4D).degreesOfFreedom();
          ++segmFromDt;
        }
        // hits from CSC segments
        if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) {
          hitsFromSegmCsc_+=(*segment)->recHits().size();
	  hitsFromSegmCsc[segmFromCsc] = hitsFromSegmCsc_;
	  chi2FromSegmCsc[segmFromCsc] = (*segment)->chi2();
	  ndofFromSegmCsc[segmFromCsc] = (*segment)->degreesOfFreedom();
          segmFromCsc++;
        }

      }    // end for muonTransientTrackingRecHit                
      
      
      // hits from track (the part concerning the rechits of rpc. info from Davide)
      //bool fbarrel[6], fendcap[3]; // flag of layers crossed
      //for (int i=0; i<6; i++) fbarrel[i] = false;
      //for (int i=0; i<3; i++) fendcap[i] = false;

      typedef std::map<RecHitIter, RecHitIter> RecToRecHitMap;
      RecToRecHitMap refToRecHitMapGLB;

      vector<RPCDetId> rpcIdGlbHit;

      for(trackingRecHit_iterator recHit = ((reco::Track)*glbOfGlobalRef).recHitsBegin(); recHit != ((reco::Track)*glbOfGlobalRef).recHitsEnd(); ++recHit){
	if (!(*recHit)->isValid()) continue;
	if ((*recHit)->isValid()) validHitsFromTrack++;
	hitsFromTrack++;
	DetId id = (*recHit)->geographicalId();

	// hits from DT
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
	  hitsFromDt++;   
	  /*
	  const DTChamber * d1;
	  d1 = dtGeometry->chamber(id);
	  DTChamberId dtId = d1->id();
	  int wheel = dtId.wheel();
	  int sector = dtId.sector();
	  int station = dtId.station(); 
	  */
	}
	// hits from CSC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) 
	  hitsFromCsc++;

	// hits from RPC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::RPC ) {
	  if(Debug_) cout<< " mark 3 "<<endl;
	  LocalPoint l;
	  l = (*recHit)->localPosition();
	  const RPCRoll * r1;
	  r1 = rpcGeometry->roll(id);
	  GlobalPoint G = r1->surface().toGlobal(l);
	  rpcGlobalPoints.push_back(G);
	  RPCDetId rpcId = r1->id();
	  rpcIdGlbHit.push_back(rpcId);

	  region[hitsFromRpc]= rpcId.region();
	  ring[hitsFromRpc] = rpcId.ring();
	  sector[hitsFromRpc] = rpcId.sector();
	  subsector[hitsFromRpc] = rpcId.subsector();
	  station[hitsFromRpc] = rpcId.station(); // 1-4 for barrel, disk 1-3 for endcaps 
	  sublayer[hitsFromRpc] = rpcId.layer();    // 1 (in) or 2 (out) just for barrel RB1 and RB2
	  recX[hitsFromRpc] = l.x();

	  hitsFromRpc++;
	}// End loop of RPC hits for GLB

	/*
	// hits from Tracker
	if (id.det() == DetId::Tracker){
	  hitsFromTk++;
	  //cout << "Tracker hits " << hitsFromTk << endl;
	}
	*/
      } // End loop of GLB muon rechits 


      TrackRef staOfGlobalRef = (*muons)[muidx].standAloneMuon();
      eta_STA = staOfGlobalRef->eta();
      phi_STA = staOfGlobalRef->phi();
      pt_STA = staOfGlobalRef->pt();
      /*
      float p_STA = staOfGlobalRef->p();
      float chi2_STA = staOfGlobalRef->chi2();
      float normChi2_STA = staOfGlobalRef->normalizedChi2();
      */
      // hit counters
      hitsFromRpcSTA=0;
      hitsFromDtSTA=0;
      hitsFromCscSTA=0;

      // hits from STA muon track
      for(trackingRecHit_iterator recHit = ((reco::Track)*staOfGlobalRef).recHitsBegin(); recHit != ((reco::Track)*staOfGlobalRef).recHitsEnd(); ++recHit){
	DetId id = (*recHit)->geographicalId();
	if (!(*recHit)->isValid()) continue;

	// hits from DT
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
	  hitsFromDtSTA++;   
	  /*
	  const DTChamber * d1;
	  d1 = dtGeometry->chamber(id);
	  DTChamberId dtId = d1->id();
	  int wheel = dtId.wheel();
	  int sector = dtId.sector();
	  int station = dtId.station(); 
	  */
	}
	// hits from CSC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) 
	  hitsFromCscSTA++;

	// hits from RPC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::RPC ) {
	  const RPCRoll * r1; LocalPoint l;
          l = (*recHit)->localPosition();
	  r1 = rpcGeometry->roll(id);
	  //GlobalPoint G = r1->surface().toGlobal(l);
	  RPCDetId rpcId = r1->id();

	  regionSTA[hitsFromRpcSTA] = rpcId.region();
	  ringSTA[hitsFromRpcSTA] = rpcId.ring();
	  sectorSTA[hitsFromRpcSTA] = rpcId.sector();
	  subsectorSTA[hitsFromRpcSTA] = rpcId.subsector();
	  stationSTA[hitsFromRpcSTA] = rpcId.station(); // 1-4 for barrel, disk 1-3 for endcaps 
	  sublayerSTA[hitsFromRpcSTA] = rpcId.layer();    // 1 (in) or 2 (out) just for barrel RB1 and RB2
	  recXSTA[hitsFromRpcSTA] = l.x();

	  hitsFromRpcSTA++;

	} // End loop of RPC hit for STA
      } // End loop of STA muon rechits

      t1->Fill();
    } // end if subset of muons



    ////////////////////////////////////////
    /////////// For tracker muons /////////
    ///////////////////////////////////////

    if (isTrackerMu[muidx] && abs(eta)<1.8) {
      TrackRef trkOfGlobalRef = (*muons)[muidx].track();
      //TrackRef trkOfGlobalRef = (*muons)[muidx].globalTrack();
      pt_track = trkOfGlobalRef->pt();
      eta_track = trkOfGlobalRef->eta();
      phi_track = trkOfGlobalRef->phi();

      // get segments from track
      MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(event, eventSetup, (reco::Track)*trkOfGlobalRef);
      segmFromDt_track=0;
      segmFromCsc_track=0;

      for (MuonTransientTrackingRecHit::MuonRecHitContainer::const_iterator segment=segments.begin(); segment!=segments.end(); segment++) {
        DetId id = (*segment)->geographicalId();
	int hitsFromSegmDt_track_=0;
	int hitsFromSegmCsc_track_=0;	

        // hits from DT segments
        if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
          const DTRecSegment4D *seg4D = dynamic_cast<const DTRecSegment4D*>((*segment)->hit());
          if((*seg4D).hasPhi())
            hitsFromSegmDt_track_+=(*seg4D).phiSegment()->specificRecHits().size();
          if((*seg4D).hasZed())
            hitsFromSegmDt_track_+=(*seg4D).zSegment()->specificRecHits().size();
	  hitsFromSegmDt_track[segmFromDt_track] = hitsFromSegmDt_track_;
	  chi2FromSegmDt_track[segmFromDt_track] = (*seg4D).chi2(); 
	  ndofFromSegmDt_track[segmFromDt_track] = (*seg4D).degreesOfFreedom();
          ++segmFromDt_track;
        }
        // hits from CSC segments
        if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) {
          hitsFromSegmCsc_track_+=(*segment)->recHits().size();
	  hitsFromSegmCsc_track[segmFromCsc_track] = hitsFromSegmCsc_track_;
	  chi2FromSegmCsc_track[segmFromCsc_track] = (*segment)->chi2();
	  ndofFromSegmCsc_track[segmFromCsc_track] = (*segment)->degreesOfFreedom();
          segmFromCsc_track++;
        }

      }    // end for muonTransientTrackingRecHit                

      // hit counters
      hitsFromRpc_track=0;
      hitsFromDt_track=0;
      hitsFromCsc_track=0;

      // hits from _track muon track
      for(trackingRecHit_iterator recHit = ((reco::Track)*trkOfGlobalRef).recHitsBegin(); recHit != ((reco::Track)*trkOfGlobalRef).recHitsEnd(); ++recHit){
	DetId id = (*recHit)->geographicalId();
	if (!(*recHit)->isValid()) continue;

	// hits from DT
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT ) {
	  hitsFromDt_track++;   
	}
	// hits from CSC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC ) 
	  hitsFromCsc_track++;

	// hits from RPC
	if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::RPC ) {
	  const RPCRoll * r1; LocalPoint l;
          l = (*recHit)->localPosition();
	  r1 = rpcGeometry->roll(id);
	  //GlobalPoint G = r1->surface().toGlobal(l);
	  RPCDetId rpcId = r1->id();

	  region_track[hitsFromRpc_track] = rpcId.region();
	  ring_track[hitsFromRpc_track] = rpcId.ring();
	  sector_track[hitsFromRpc_track] = rpcId.sector();
	  subsector_track[hitsFromRpc_track] = rpcId.subsector();
	  station_track[hitsFromRpc_track] = rpcId.station(); // 1-4 for barrel, disk 1-3 for endcaps 
	  sublayer_track[hitsFromRpc_track] = rpcId.layer();    // 1 (in) or 2 (out) just for barrel RB1 and RB2
	  recX_track[hitsFromRpc_track] = l.x();

	  hitsFromRpc_track++;

	} // End loop of RPC hit for Tracker Muons
      } // End loop of Tracker Muon rechits
      t3->Fill();
    } // end if subset of muon==isTrackerMuon


  }  // end loop on muons

  t0->Fill();


  // loop over all CSC segments                                                                                                                         
  int icscsg=0;
  for (CSCSegmentCollection::const_iterator segmentCSC = allCSCSegments->begin();
       segmentCSC!=allCSCSegments->end(); ++segmentCSC){
    icscsg++;
    CSCDetId id  = (CSCDetId)(*segmentCSC).cscDetId();
    endcap  = id.endcap();
    ringCSC    = id.ring();
    stationCSC = id.station();
    chamberCSC = id.chamber();
    // Get pointer to the chamber:                                                                                                                    
    const CSCChamber* cscchamber = cscGeometry->chamber(id);
    
    LocalPoint CSCLocalPosition = segmentCSC->localPosition();
    GlobalPoint CSCGlobalPosition = cscchamber->toGlobal(CSCLocalPosition);

    globalX = CSCGlobalPosition.x();
    globalY = CSCGlobalPosition.y();
    //cout << "CSC XY: "<<CSCGlobalPosition.x()<<", "<< CSCGlobalPosition.y() << endl;
    t2->Fill();
  }
  
  
 }


// ------------ method called once each job just before starting event loop  ------------
void RecHitAnal::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void RecHitAnal::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnal);
