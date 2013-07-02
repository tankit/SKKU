#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;
using namespace edm;
using namespace reco;

class PromptMuonSelector : public edm::EDFilter
{
public:
  PromptMuonSelector(const edm::ParameterSet& pset);
  virtual ~PromptMuonSelector() {};
  void beginJob() {};
  void endJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag pvLabel_;
  bool doBeamSpot_;
  edm::InputTag muonLabel_;
  double maxDxy_;
  double maxDz_;

};

PromptMuonSelector::PromptMuonSelector(const edm::ParameterSet& pset)
{
  if ( pset.exists("vertex") )
  {
    pvLabel_ = pset.getParameter<edm::InputTag>("vertex");
    doBeamSpot_ = false;
  }
  else
  {
    pvLabel_ = pset.getParameter<edm::InputTag>("beamSpot");
    doBeamSpot_ = true;
  }

  muonLabel_ = pset.getParameter<edm::InputTag>("src");
  maxDxy_ = pset.getUntrackedParameter<double>("maxDxy", 0.2);
  maxDz_ = pset.getUntrackedParameter<double>("maxDz", 0.5);

  produces<std::vector<reco::Muon> >("");
}

bool PromptMuonSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  math::XYZPoint pvPos;
  if ( doBeamSpot_ )
  {
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    event.getByLabel(pvLabel_, beamSpotHandle);
    pvPos = beamSpotHandle->position();
  }
  else
  {
    edm::Handle<reco::VertexCollection> pvHandle;
    event.getByLabel(pvLabel_, pvHandle);
    pvPos = pvHandle->at(0).position();
  }

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  if ( !muonHandle.isValid() ) return false;

  std::auto_ptr<std::vector<reco::Muon> > promptMuons(new std::vector<reco::Muon>);

  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {
    bool isGoodIP = true;
    if ( muon->innerTrack().isAvailable() ) 
    {
      if ( muon->innerTrack()->dxy(pvPos) > maxDxy_ ) isGoodIP = false;
      if ( abs(muon->innerTrack()->dz(pvPos)) > maxDz_ ) isGoodIP = false;
    }
    else
    {
      if ( not muon->outerTrack().isAvailable() ) 
      {
        LogVerbatim("PromptMuonSelector") << "No referenced tracks are available, skip the muon";
        continue;
      }
      if ( muon->outerTrack()->dxy(pvPos) > maxDxy_ ) isGoodIP = false;
      if ( abs(muon->outerTrack()->dz(pvPos)) > maxDz_ ) isGoodIP = false;
    }

    if ( isGoodIP )
    {
      promptMuons->push_back(*muon);
    }
  }
  const bool isEmpty = promptMuons->empty();

  event.put(promptMuons);

  return !isEmpty;
}

DEFINE_FWK_MODULE(PromptMuonSelector);
