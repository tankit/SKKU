import FWCore.ParameterSet.Config as cms

### HLT filter
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
ZMuHLTFilter = copy.deepcopy(hltHighLevel)
#ZMuHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT") #default for Data
ZMuHLTFilter.throw = cms.bool(False)
#ZMuHLTFilter.HLTPaths = ["HLT_Mu*","HLT_IsoMu*","HLT_DoubleMu*"]
ZMuHLTFilter.HLTPaths = ["HLT_Mu*","HLT_IsoMu*"]

### Z -> MuMu candidates

# Get muons of needed quality for Zs
looseMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("muons"),
                             #cut = cms.string('pt > 10 && abs(eta)<2.4 && isGlobalMuon = 1 && isTrackerMuon = 1 && abs(innerTrack().dxy)<2.0'),
                             cut = cms.string('pt > 10 && abs(eta)<2.4 && isGlobalMuon = 1 && abs(innerTrack().dxy)<2.0'),
                             filter = cms.bool(True)
                             )

tightMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("looseMuonsForZ"),
                             cut = cms.string('pt > 15 && globalTrack().normalizedChi2<10.0 && isolationR03().sumPt<3.0 && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt)<0.2*pt && track().hitPattern().numberOfValidTrackerHits()>10 && numberOfMatches>1'),
                             filter = cms.bool(True)                                
                             )

generalTracksFilter = cms.EDFilter("TrackCountFilter",
                                   src = cms.InputTag('generalTracks'),
                                   #cut = cms.string('pt > 5 && abs(eta)<2.4'),
                                   cut = cms.string('pt > 15 && abs(eta)<2.4'),
                                   minNumber = cms.uint32(2) 
                                   )

goodTracks = cms.EDFilter("TrackSelector",
                          src = cms.InputTag("generalTracks"), # or cms.InputTag("standAloneMuons","UpdatedAtVtx"),
                          cut = cms.string(""),
                          #cut = cms.string("numberOfValidHits >= 10 && normalizedChi2 < 5 && abs(d0) < 2 && abs(dz) < 30"),
                          #cut = cms.string("numberOfValidHits >= 10"),
                          )

trackCands  = cms.EDProducer("ConcreteChargedCandidateProducer",
                             src  = cms.InputTag("goodTracks"),
                             particleType = cms.string("mu+"),     # this is needed to define a mass
                             )

trackProbes = cms.EDFilter("CandViewRefSelector",
                           src = cms.InputTag("trackCands"),
                           cut = cms.string("pt>15 && abs(eta)<1.8"),
                           #cut = cms.string("pt>5 && abs(eta)<1.8"),
                           )

# build Z-> MuMu candidates
dimuons = cms.EDProducer("CandViewShallowCloneCombiner",
                         checkCharge = cms.bool(False),
                         cut = cms.string('mass > 30'),
                         decay = cms.string("tightMuonsForZ@+ trackProbes@-") # charge conjugate states are implied; 'tagMuons' and 'trkProbes' should be collections of Candidates
                         )

# Z filter
dimuonsFilter = cms.EDFilter("CandViewCountFilter",
                             src = cms.InputTag("dimuons"),
                             minNumber = cms.uint32(1)
                             )

# TagMu Skim sequence
TagMuSelSeq = cms.Sequence(ZMuHLTFilter *
                           looseMuonsForZ *
                           tightMuonsForZ * #generalTracksFilter *
                           goodTracks *
                           trackCands *
                           trackProbes *
                           dimuons *
                           dimuonsFilter
                           )


