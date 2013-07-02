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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
  edm::InputTag pvLabel_;
  bool doBeamSpot_;
  edm::InputTag trackCandLabel_;
  double maxDxy_, maxDz_;

};

PromptTrackCandSelector::PromptTrackCandSelector(const edm::ParameterSet& pset)
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

  trackCandLabel_ = pset.getParameter<edm::InputTag>("src");
  maxDxy_ = pset.getUntrackedParameter<double>("maxDxy", 0.2);
  maxDz_ = pset.getUntrackedParameter<double>("maxDz", 0.2);

  produces<std::vector<reco::RecoChargedCandidate> >("");
}

bool PromptTrackCandSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
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

  edm::Handle<edm::View<reco::RecoChargedCandidate> > trackCandHandle;
  event.getByLabel(trackCandLabel_, trackCandHandle);
  if ( !trackCandHandle.isValid() ) return false;

  std::auto_ptr<std::vector<reco::RecoChargedCandidate> > promptTrackCands(new std::vector<reco::RecoChargedCandidate>);

  for ( edm::View<reco::RecoChargedCandidate>::const_iterator trackCand = trackCandHandle->begin();
        trackCand != trackCandHandle->end(); ++trackCand )
  {
    if ( trackCand->track()->dxy(pvPos) > maxDxy_ ) continue;
    if ( abs(trackCand->track()->dz(pvPos)) > maxDz_ ) continue;

    promptTrackCands->push_back(*trackCand);
  }

  const bool isEmpty = promptTrackCands->empty();

  event.put(promptTrackCands);

  return !isEmpty;
}

DEFINE_FWK_MODULE(PromptTrackCandSelector);
