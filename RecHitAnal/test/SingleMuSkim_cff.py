import FWCore.ParameterSet.Config as cms

### HLT filter
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
ZMuHLTFilter = copy.deepcopy(hltHighLevel)
ZMuHLTFilter.throw = cms.bool(False)
ZMuHLTFilter.HLTPaths = ["HLT_Mu*","HLT_IsoMu*","HLT_DoubleMu*"]

### Z -> MuMu candidates

# Get muons of needed quality for Zs
looseMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("muons"),
                             cut = cms.string('pt > 10 && abs(eta)<2.4 && isGlobalMuon = 1 && isTrackerMuon = 1 && abs(innerTrack().dxy)<2.0'),
                             filter = cms.bool(True)                                
                             )

tightMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("looseMuonsForZ"),
                             cut = cms.string('pt > 15 && globalTrack().normalizedChi2<10.0 && isolationR03().sumPt<3.0 && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt)<0.2*pt && globalTrack().hitPattern().numberOfValidTrackerHits>10'),
                             filter = cms.bool(True)                                
                             )

# build Z-> MuMu candidates
dimuons = cms.EDProducer("CandViewShallowCloneCombiner",
                         checkCharge = cms.bool(False),
                         cut = cms.string('mass > 30'),
                         decay = cms.string("tightMuonsForZ looseMuonsForZ")
                         )

# Z filter
dimuonsFilter = cms.EDFilter("CandViewCountFilter",
                             src = cms.InputTag("dimuons"),
                             minNumber = cms.uint32(1)
                             )

generalTracksFilter = cms.EDFilter("TrackCountFilter",
                                   src = cms.InputTag('generalTracks'),
                                   cut = cms.string('pt > 5'),
                                   minNumber = cms.uint32(2)
                                   )

# Z Skim sequence
diMuonSelSeq = cms.Sequence(ZMuHLTFilter *
                            looseMuonsForZ *
                            tightMuonsForZ *
                            generalTracksFilter
                            )


