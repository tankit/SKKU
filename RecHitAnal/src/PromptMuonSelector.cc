#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
  edm::InputTag beamSpotLabel_;
  edm::InputTag muonLabel_;
  double maxDxy_;

};

PromptMuonSelector::PromptMuonSelector(const edm::ParameterSet& pset)
{
  beamSpotLabel_ = pset.getParameter<edm::InputTag>("beamSpot");
  muonLabel_ = pset.getParameter<edm::InputTag>("muon");
  maxDxy_ = pset.getUntrackedParameter<double>("maxDxy", 0.2);

  produces<std::vector<reco::Muon> >("");
}

bool PromptMuonSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel(beamSpotLabel_, beamSpotHandle);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  if ( !beamSpotHandle.isValid() ) return false;
  if ( !muonHandle.isValid() ) return false;

  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  std::auto_ptr<std::vector<reco::Muon> > promptMuons(new std::vector<reco::Muon>);

  for ( edm::View<reco::Muon>::const_iterator muon = muonHandle->begin();
        muon != muonHandle->end(); ++muon )
  {
    if ( abs(muon->innerTrack()->dxy(beamSpot)) < maxDxy_ )
    {
      promptMuons->push_back(*muon);
    } 
  }

  event.put(promptMuons);

  if ( promptMuons->empty() ) return false;
  return true;
}

DEFINE_FWK_MODULE(PromptMuonSelector);
