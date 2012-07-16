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
  double maxDz_;

};

PromptMuonSelector::PromptMuonSelector(const edm::ParameterSet& pset)
{
  beamSpotLabel_ = pset.getParameter<edm::InputTag>("beamSpot");
  muonLabel_ = pset.getParameter<edm::InputTag>("src");
  maxDxy_ = pset.getUntrackedParameter<double>("maxDxy", 0.2);
  maxDz_ = pset.getUntrackedParameter<double>("maxDz", 0.5);

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
    if ( not muon->innerTrack().isAvailable() ) {
      LogVerbatim("PromptMuonSelector") << "No inner tracks are available, skim the muon";
      continue;
    }
    if ( abs(muon->innerTrack()->dxy(beamSpot)) < maxDxy_ )
    //if ( abs(muon->innerTrack()->dxy(beamSpot)) < maxDxy_ && abs(muon->innerTrack()->dz(beamSpot)) < maxDz_ )
    {
      promptMuons->push_back(*muon);
    } 
  }
  const bool isEmpty = promptMuons->empty();

  event.put(promptMuons);

  return !isEmpty;
}

DEFINE_FWK_MODULE(PromptMuonSelector);
