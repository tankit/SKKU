#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

//#include "DataFormats/CandidateReco/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;
using namespace edm;
using namespace reco;

class PromptTrackCandSelector : public edm::EDFilter
{
public:
  PromptTrackCandSelector(const edm::ParameterSet& pset);
  virtual ~PromptTrackCandSelector() {};
  void beginJob() {};
  void endJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag beamSpotLabel_;
  edm::InputTag trackCandLabel_;
  double maxDxy_;

};

PromptTrackCandSelector::PromptTrackCandSelector(const edm::ParameterSet& pset)
{
  beamSpotLabel_ = pset.getParameter<edm::InputTag>("beamSpot");
  trackCandLabel_ = pset.getParameter<edm::InputTag>("src");
  maxDxy_ = pset.getUntrackedParameter<double>("maxDxy", 0.2);

  produces<std::vector<reco::RecoChargedCandidate> >("");
}

bool PromptTrackCandSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel(beamSpotLabel_, beamSpotHandle);

  edm::Handle<edm::View<reco::RecoChargedCandidate> > trackCandHandle;
  event.getByLabel(trackCandLabel_, trackCandHandle);

  if ( !beamSpotHandle.isValid() ) return false;
  if ( !trackCandHandle.isValid() ) return false;

  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  std::auto_ptr<std::vector<reco::RecoChargedCandidate> > promptTrackCands(new std::vector<reco::RecoChargedCandidate>);

  for ( edm::View<reco::RecoChargedCandidate>::const_iterator trackCand = trackCandHandle->begin();
        trackCand != trackCandHandle->end(); ++trackCand )
  {
    if ( abs(trackCand->track()->dxy(beamSpot)) < maxDxy_ )
    {
      promptTrackCands->push_back(*trackCand);
    } 
  }

  const bool isEmpty = promptTrackCands->empty();

  event.put(promptTrackCands);

  return !isEmpty;
}

DEFINE_FWK_MODULE(PromptTrackCandSelector);
