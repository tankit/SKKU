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
// $Id: RecHitAnal.cc,v 1.3 2011/07/06 15:40:19 hkseo Exp $=======
// $Id: RecHitAnal.cc,v 1.4 2011/07/27 14:26:23 hkseo Exp $>>>>>>> 1.4
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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
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
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DQM/RPCMonitorDigi/interface/utils.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
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
#include<iostream>
#include<cmath>
#include<assert.h>
#include "TPaletteAxis.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TFolder.h"
#include "TTree.h"
#include "TFile.h"
#include"TH1.h"
#include "TH2F.h"
#include "TLorentzVector.h"


using namespace edm; 
using namespace std;
using namespace reco;

static const int kMax = 10;

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

  TTree *t0; //Tree for each event
  Int_t Run;
  Int_t eventNumber;
  Int_t totalEvents;
  Float_t angle;
  double d0;
  Int_t nTracks;
  Int_t nMuons;
  
  TTree *t1; // Tree for each muon
  Int_t hitsFromRpc, hitsFromDt, hitsFromCsc;
  Int_t hitsFromRpcSTA, hitsFromDtSTA, hitsFromCscSTA;
  Float_t eta, phi, pt;
  Float_t eta_STA, phi_STA, pt_STA;

  //for GLB
  Int_t region[kMax];
  Int_t ring[kMax];
  Int_t sector[kMax];
  Int_t subsector[kMax];
  Int_t station[kMax];  // 1-4 for barrel, disk 1-3 for endcaps 
  Int_t sublayer[kMax]; // 1 (in) or 2 (out) just for barrel RB1 and RB2  
  Float_t recX[kMax];
  //for STA
  Int_t regionSTA[kMax];
  Int_t ringSTA[kMax];
  Int_t sectorSTA[kMax];
  Int_t subsectorSTA[kMax];
  Int_t stationSTA[kMax];  // 1-4 for barrel, disk 1-3 for endcaps 
  Int_t sublayerSTA[kMax]; // 1 (in) or 2 (out) just for barrel RB1 and RB2  
  Float_t recXSTA[kMax];

  bool Debug_;

  InputTag RPCRecHits_;
  InputTag DT4DSegments_;
  InputTag CSCSegments_;
  InputTag muons_;
  InputTag StandAloneMuons_;
  InputTag genTracksLabel_;

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

  //const ParameterSet SegmentsTrackAssociatorParameters = cfg.getParameter<ParameterSet>("SegmentsTrackAssociatorParameters");
  //theSegmentsAssociator = new SegmentsTrackAssociator(SegmentsTrackAssociatorParameters);
  
  //debug
  Debug_ = cfg.getUntrackedParameter<bool>("Debug",false);
  
  //RPC hits and digis 
  RPCRecHits_ = cfg.getParameter<InputTag>("RPCRecHits"); 
  DT4DSegments_ = cfg.getParameter<InputTag>("DT4DSegments");
  CSCSegments_ = cfg.getParameter<InputTag>("CSCSegments");
  muons_ = cfg.getParameter<InputTag>("muons");
  StandAloneMuons_ = cfg.getParameter<InputTag>("StandAloneMuons");
  genTracksLabel_ = cfg.getParameter<InputTag>("generalTracks");
  RPCRecHits_ = cfg.getParameter<InputTag>("RPCRecHits"); 


  //create the Tree and branches
  t0 = fs->make<TTree>("Events","");
  t0->Branch("Run",         &Run,              "Run/I");
  t0->Branch("eventNumber", &eventNumber,      "eventNumber/I");
  t0->Branch("totalEvents", &totalEvents,      "totalEvents/I");
  t0->Branch("angle",       &angle,            "angle/F");
  t0->Branch("d0",          &d0,               "d0/D");
  t0->Branch("nTracks",     &nTracks,          "nTracks/I");
  t0->Branch("nMuons",      &nMuons,           "nMuons/I");

  t1 = fs->make<TTree>("muon","");
  t1->Branch("Run",         &Run,              "Run/I");
  t1->Branch("eventNumber", &eventNumber,      "eventNumber/I");
  t1->Branch("angle",       &angle,            "angle/F");
  t1->Branch("d0",          &d0,               "d0/D");
  t1->Branch("nTracks",     &nTracks,          "nTracks/I");
  t1->Branch("nMuons",      &nMuons,           "nMuons/I");
  t1->Branch("nRpcHit",     &hitsFromRpc,      "nRpcHit/I");
  t1->Branch("nDtHit",      &hitsFromDt,       "nDtHit/I");
  t1->Branch("nCscHit",     &hitsFromCsc,      "nCscHit/I");
  t1->Branch("eta",         &eta,              "eta/F");
  t1->Branch("phi",         &phi,              "phi/F");
  t1->Branch("pt",          &pt,               "pt/F");
  t1->Branch("nRpcHitSTA",  &hitsFromRpcSTA,   "nRpcHitSTA/I");
  t1->Branch("nDtHitSTA",   &hitsFromDtSTA,    "nDtHitSTA/I");
  t1->Branch("nCscHitSTA",  &hitsFromCscSTA,   "nCscHitSTA/I");
  t1->Branch("etaSTA",      &eta_STA,          "eta/F");
  t1->Branch("phiSTA",      &phi_STA,          "phi/F");
  t1->Branch("ptSTA",       &pt_STA,           "pt/F");

  // For global muons
  t1->Branch("region",      region,            "region[nRpcHit]/I");
  t1->Branch("ring",        ring,              "ring[nRpcHit]/I");
  t1->Branch("sector",      sector,            "sector[nRpcHit]/I");
  t1->Branch("subsector",   subsector,         "subsector[nRpcHit]/I");
  t1->Branch("station",     station,           "station[nRpcHit]/I");
  t1->Branch("sublayer",    sublayer,          "sublayer[nRpcHit]/I");
  t1->Branch("localX",      recX,              "localX[nRpcHit]/F");
  // For standalone muons
  t1->Branch("regionSTA",   regionSTA,         "regionSTA[nRpcHitSTA]/I");
  t1->Branch("ringSTA",     ringSTA,           "ringSTA[nRpcHitSTA]/I");
  t1->Branch("sectorSTA",   sectorSTA,         "sectorSTA[nRpcHitSTA]/I");
  t1->Branch("subsectorSTA",subsectorSTA,      "subsectorSTA[nRpcHitSTA]/I");
  t1->Branch("stationSTA",  stationSTA,        "stationSTA[nRpcHitSTA]/I");
  t1->Branch("sublayerSTA", sublayerSTA,       "sublayerSTA[nRpcHitSTA]/I");
  t1->Branch("localXSTA",   recXSTA,           "localXSTA[nRpcHitSTA]/F");
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

  int nCSCseg = allCSCSegments->size();
  int nDTseg = allDT4DSegments->size();
  int nRPCRecHits = rpcRecHits->size();
  nMuons = muons->size();
  int nStaMuons = StandAloneMuons->size();
  nTracks = generalTracks->size();

  if (nMuons<1) return;   // for now ask for at least one muon


  //*****************************************************************
  //*  Calculate 3D angle for muons in order to remove cosmic muons *
  //*****************************************************************

  bool isTrackerMu = false;
  bool isStaMu = false;
  bool isGlobalMu = false;
  vector<TLorentzVector> zmumu;
  angle = -1;

  for (int muidx=0; muidx< nMuons; muidx++) {  // loop on all muons
    isTrackerMu = false;
    isStaMu = false;
    isGlobalMu = false;
    float eta_ = ((*muons)[muidx]).eta();
    //float p_ = ((*muons)[muidx]).p();
    float px_ = ((*muons)[muidx]).px();
    float py_ = ((*muons)[muidx]).py();
    float pz_ = ((*muons)[muidx]).pz();
    float energy_ = ((*muons)[muidx]).energy();
    //int charge_ = ((*muons)[muidx]).charge();
    
    if ((*muons)[muidx].isTrackerMuon() && abs(eta_)<2.4) isTrackerMu = true;
    if ((*muons)[muidx].isStandAloneMuon() && abs(eta_)<2.4) isStaMu = true;
    if ((*muons)[muidx].isGlobalMuon() && abs(eta_)<2.4) isGlobalMu = true;

    if (isGlobalMu) zmumu.push_back( TLorentzVector(px_,py_,pz_,energy_) );
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
    isTrackerMu = false;
    isStaMu = false;
    isGlobalMu = false;
    // vector of RpcRecHits associated to sta muon
    //    vector<RPCRecHit> RPCHistOnSta_List;
    
    eta = ((*muons)[muidx]).eta();
    phi = ((*muons)[muidx]).phi();
    pt = ((*muons)[muidx]).pt();
    //float p = ((*muons)[muidx]).p();
    //float charge = ((*muons)[muidx]).charge();

    if ((*muons)[muidx].isTrackerMuon() && abs(eta)<1.6) isTrackerMu = true;
    if ((*muons)[muidx].isStandAloneMuon() && abs(eta)<1.6) isStaMu = true;
    if ((*muons)[muidx].isGlobalMuon() && abs(eta)<1.6) isGlobalMu = true;


    if (isGlobalMu && isStaMu && abs(eta)<1.6) {
      if(Debug_) cout<< " mark 2 "<<endl;
      //initialize
      d0 = 0;

      TrackRef glbOfGlobalRef = (*muons)[muidx].globalTrack();
      //TrackRef trkOfGlobalRef = (*muons)[muidx].track();

      /*
      float pt_track = trkOfGlobalRef->pt();
      float p_track = trkOfGlobalRef->p();
      float eta_track = trkOfGlobalRef->eta();
      float phi_track = trkOfGlobalRef->phi();
      float chi2_track = trkOfGlobalRef->chi2();
      float normChi2_track = trkOfGlobalRef->normalizedChi2();
      */

      math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
      d0 = -1.* glbOfGlobalRef->dxy(point);
      //cout << "dxy: "<<glbOfGlobalRef->dxy()<<" "<<d0<< endl;

      // temporary part
      //      pt = pt_track;
      //      p = p_track;
      //      eta = eta_track;

      // get segments from track
      //MuonTransientTrackingRecHit::MuonRecHitContainer segments = theSegmentsAssociator->associate(event, eventSetup, (reco::Track)*staOfGlobalRef );

      // hit counters
      hitsFromDt=0;
      hitsFromCsc=0;
      hitsFromRpc=0;
      int hitsFromTk=0;
      int hitsFromTrack=0;
      int validHitsFromTrack=0;
      int hitsFromSegmDt=0;
      int hitsFromSegmCsc=0;
      // segment counters
      int segmFromDt=0;
      int segmFromCsc=0;

      // hits from track (the part concerning the rechits of rpc. info from Davide)
      bool fbarrel[6], fendcap[3]; // flag of layers crossed
      for (int i=0; i<6; i++) fbarrel[i] = false;
      for (int i=0; i<3; i++) fendcap[i] = false;

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
	  GlobalPoint G = r1->surface().toGlobal(l);
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


  }  // end loop on muons

  t0->Fill();
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
