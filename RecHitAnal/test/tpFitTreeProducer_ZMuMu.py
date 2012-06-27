import FWCore.ParameterSet.Config as cms
##                      _              _
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##
################################################
maxDxy = 0.2
MC_flag = False
#MC_flag = True
GLOBAL_TAG = 'GR_R_52_V8::All'
if MC_flag:
    GLOBAL_TAG = 'START52_V5::All'

#HLTPath1 = "HLT_IsoMu24_v1"
#HLTPath2 = "HLT_IsoMu24_v8"
#HLTPath3 = "HLT_IsoMu24_v9"
HLTProcessName = "HLT"

RECOProcess = "RECO"
JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.6 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99 && pt>15.0"

process = cms.Process("TagProbe")

######### EXAMPLE CFG 
###  A simple test of runnning T&P on Zmumu to determine muon isolation and identification efficiencies
###  More a showcase of the tool than an actual physics example

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.GlobalTag.globaltag = GLOBAL_TAG
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
         'file:/afs/cern.ch/user/m/mskim/workspace/public/rpc/00AD4245-A5B5-E111-A1E8-001EC9D8B54A.root',
         'file:/afs/cern.ch/user/m/mskim/workspace/public/rpc/0215F1C2-BBAE-E111-A01F-485B39800B69.root',
    )
)
#for line in open('Run2012B-ZMu-PromptSkim-v1_255.txt').readlines():
#    line = line.strip("'\", \n")
#    if '.root' not in line: continue
#
#    process.source.fileNames.append(line)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# Filter some runs
#process.load("RPCAnalysis.RunFilter.runfilter_cfi")
#process.runfilter.iRun = cms.untracked.int32(170249)
#process.runfilter.fRun = cms.untracked.int32(170527)
#process.runfilter.veto = cms.untracked.bool(True)

### HLT filter
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZMuHLTFilter = copy.deepcopy(hltHighLevel)
process.ZMuHLTFilter.throw = cms.bool(False)
#process.ZMuHLTFilter.HLTPaths = [HLTPath1,HLTPath2,HLTPath3]

process.HLTFilter = cms.Sequence(process.ZMuHLTFilter)


###############################################
# Re-make muon reco without RPC hits
# standalone muons
process.load("RecoMuon.StandAloneMuonProducer.standAloneMuons_cff")
process.standAloneMuonsNoRPC = process.standAloneMuons.clone()
process.standAloneMuonsNoRPC.STATrajBuilderParameters.FilterParameters.EnableRPCMeasurement = cms.bool(False)
process.standAloneMuonsNoRPC.STATrajBuilderParameters.BWFilterParameters.EnableRPCMeasurement = cms.bool(False)
# global muons
process.load("RecoMuon.GlobalMuonProducer.globalMuons_cff")
process.globalMuonsNoRPC = process.globalMuons.clone()
process.globalMuonsNoRPC.MuonCollectionLabel = cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
process.muontrackingNoRPC = cms.Sequence(
    process.standAloneMuonsNoRPC *
    process.globalMuonsNoRPC
    )

# muons
#process.load("RecoMuon.MuonIdentification.muons_cfi")
#process.muonsNoRPC = process.muons.clone()
#process.muonsNoRPC.inputCollectionLabels = cms.VInputTag(
#    cms.InputTag("generalTracks"),
#    cms.InputTag("globalMuonsNoRPC"),
#    cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
#    )
#process.muonsNoRPC.fillGlobalTrackQuality = cms.bool(False)
#
#process.muonIdProducerSequenceNoRPC = cms.Sequence(
#    process.muonsNoRPC
#    )
###############################################            

# Muon Id producer
process.load("RecoMuon.MuonIdentification.muons1stStep_cfi")
process.muonsNoRPC1stStep = process.muons1stStep.clone()
#process.muonsNoRPC1stStep.inputCollectionLabels = ['ctfWithMaterialTracksP5', 'globalMuonsNoRPC', 'standAloneMuonsNoRPC']
process.muonsNoRPC1stStep.inputCollectionLabels = cms.VInputTag(
    cms.InputTag("generalTracks"),
    cms.InputTag("globalMuonsNoRPC"),
    cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
    )
process.muonsNoRPC1stStep.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
#process.muonsNoRPC1stStep.fillIsolation = cms.bool(True) #True by default
process.muonsNoRPC1stStep.fillGlobalTrackQuality = cms.bool(False)
##process.muonsNoRPC1stStep.TrackExtractorPSet.inputTrackCollection = cms.InputTag("ctfWithMaterialTracksP5") #'ctfWithMaterialTracksP5'
##process.muonsNoRPC1stStep.CaloExtractorPSet.CenterConeOnCalIntersection = cms.bool(True)
#process.muonsNoRPC1stStep.fillGlobalTrackRefits = cms.bool(False) #True by default

process.load("RecoMuon.MuonIdentification.muons_cfi")
process.muonsNoRPC = process.muons.clone()
process.muonsNoRPC.InputMuons = cms.InputTag("muonsNoRPC1stStep")
process.muonsNoRPC.FillPFMomentumAndAssociation = False
process.muonsNoRPC.FillPFIsolation = False
process.muonsNoRPC.FillSelectorMaps = False
process.muonsNoRPC.FillCosmicsIdMap =  False
process.muonsNoRPC.FillTimingInfo = False
process.muonsNoRPC.FillDetectorBasedIsolation = False
process.muonsNoRPC.FillShoweringInfo = False

process.muonIdProducerSequenceNoRPC = cms.Sequence(
    process.muonsNoRPC1stStep *
    process.muonsNoRPC
    )


##############
# TAG MUON
##############

# Trigger  ##################
process.PassingHLT = cms.EDProducer("trgMatchedMuonProducer",
                                    InputProducer = cms.InputTag("promptMuons"),
                                    hltTags = cms.VInputTag(
                                    #cms.InputTag(HLTPath1,"", HLTProcessName),
                                    #cms.InputTag(HLTPath2,"", HLTProcessName),
                                    #cms.InputTag(HLTPath3,"", HLTProcessName)
                                    ),
                                    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
                                    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
                                    )

for version in range(1,21):
    process.PassingHLT.hltTags.append(cms.InputTag("HLT_IsoMu24_v%d::%s" % (version, HLTProcessName)))
    process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu24_v%d" % version)

## Tags. In a real analysis we should require that the tag muon fires the trigger, 
##       that's easy with PAT muons but not RECO/AOD ones, so we won't do it here
##       (the J/Psi example shows it)
#PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % ("HLT_Mu30_v7",);

process.promptMuons = cms.EDFilter("PromptMuonSelector",
    src = cms.InputTag("muons"),
    maxDxy = cms.untracked.double(maxDxy),
    beamSpot = cms.InputTag("offlineBeamSpot"),
)
process.promptMuonsNoRPC = process.promptMuons.clone(src = cms.InputTag("muonsNoRPC"))

process.tightMuons = cms.EDFilter("MuonRefSelector",
                                src = cms.InputTag("promptMuons"),
                                cut = cms.string("isGlobalMuon && isTrackerMuon && isolationR03().sumPt<3.0"
                                                 "&& pt > 20 && abs(eta) < 2.4"
                                                 "&& track().hitPattern().numberOfValidPixelHits() > 0" #--added by Minsuk on Feb 17, 2012
                                                 "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                 "&& innerTrack().numberOfValidHits()>10 && globalTrack().normalizedChi2()<10.0"
                                                 "&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                 "&& numberOfMatchedStations>1" #--updated from numberOfMatches by Minsuk on May 12, 2012
                                                 "&& (isolationR03().sumPt+isolationR03().emEt+isolationR03().hadEt)<0.1*pt"
                                                 ), 
                                )

process.tagMuons = process.PassingHLT.clone()
process.tagMuons.InputProducer = cms.InputTag("tightMuons")

## Probes. Now we just use Tracker Muons as probes
process.probeMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("promptMuons"),
    #cut = cms.string("isTrackerMuon && pt > 10"), 
    #cut = cms.string("isTrackerMuon && pt > 20 && abs(eta)<1.8 && innerTrack.numberOfValidHits() >= 10"), 
    ##cut = cms.string("isRPCMuon && pt > 20 && abs(eta)<1.8 && innerTrack.numberOfValidHits() >= 10"),
    cut = cms.string("isGlobalMuon && pt > 20 && abs(eta) < 1.6"
                     ##"&& isolationR03().sumPt<3.0 && (isolationR03().sumPt+isolationR03().emEt+isolationR03().hadEt)<0.1*pt"
                     #"&& numberOfMatchedStations > 1"
                     #"&& pt > 0 && abs(eta) < 2.4"
                     "&& globalTrack().normalizedChi2()<10.0"
                     #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                     "&& track().hitPattern().numberOfValidPixelHits() > 0"
                     "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                     "&& innerTrack().numberOfValidHits()>10"
                     "&& abs(innerTrack().dxy)<2.0" #--added because of no IP cut from Jun 12, 2012
                     ),
)


#######################
# PROBE GENERAL TRACK
#######################

process.goodTracks = cms.EDFilter("TrackSelector",
                                  src = cms.InputTag("generalTracks"), # or cms.InputTag("standAloneMuons","UpdatedAtVtx"),
                                  #cut = cms.string("numberOfValidHits >= 10 && normalizedChi2 < 5 && abs(d0) < 2 && abs(dz) < 30"),
                                  cut = cms.string("pt>20 && abs(eta)<2.4 && numberOfValidHits >= 10"),
                                  )

process.trackCands  = cms.EDProducer("ConcreteChargedCandidateProducer",
                                     src  = cms.InputTag("goodTracks"),
                                     particleType = cms.string("mu+"),     # this is needed to define a mass
                                     )

process.promptTrackCands = cms.EDFilter("PromptTrackCandSelector",
    src = cms.InputTag("trackCands"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    maxDxy = cms.untracked.double(maxDxy),
)

process.trackProbes = cms.EDFilter("CandViewRefSelector",
                                   src = cms.InputTag("promptTrackCands"),
                                   cut = cms.string("abs(eta)<1.8"),
                                   )


############
# MUON ID
############

process.MediumMuons = cms.EDFilter("MuonRefSelector",
                                  src = cms.InputTag("promptMuons"),
                                  cut = cms.string("isGlobalMuon"
                                                   #"&& numberOfMatchedStations > 1"
                                                   "&& pt > 0 && abs(eta) < 2.4"
                                                   "&& globalTrack().normalizedChi2()<10.0"
                                                   #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                   "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                   "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                   "&& innerTrack().numberOfValidHits()>10"
                                                   ), 
                                  )


process.MediumMuonsNoRPC = cms.EDFilter("MuonRefSelector",
                                       src = cms.InputTag("promptMuonsNoRPC"),
                                       cut = cms.string("isGlobalMuon"
                                                        #"&& numberOfMatchedStations > 1"
                                                        "&& pt > 0 && abs(eta) < 2.4"
                                                        "&& globalTrack().normalizedChi2()<10.0"
                                                        #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                        "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                        "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                        "&& innerTrack().numberOfValidHits()>10"
                                                        ), 
                                       )


process.TightMuons = cms.EDFilter("MuonRefSelector",
                                  src = cms.InputTag("promptMuons"),
                                  cut = cms.string("isGlobalMuon"
                                                   "&& numberOfMatchedStations > 1"
                                                   "&& pt > 0 && abs(eta) < 2.4"
                                                   "&& globalTrack().normalizedChi2()<10.0"
                                                   "&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                   "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                   "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                   "&& innerTrack().numberOfValidHits()>10"
                                                   ),
                                   )


process.TightMuonsNoRPC = cms.EDFilter("MuonRefSelector",
                                       src = cms.InputTag("promptMuonsNoRPC"),
                                       cut = cms.string("isGlobalMuon"
                                                        "&& numberOfMatchedStations > 1"
                                                        "&& pt > 0 && abs(eta) < 2.4"
                                                        "&& globalTrack().normalizedChi2()<10.0"
                                                        "&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                        "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                        "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                        "&& innerTrack().numberOfValidHits()>10"
                                                        ),
                                        )



########################
# MATCH TRACK AND MUONS
########################

process.tkToMediumMuons = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("promptTrackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("MediumMuons"), # to all global muons
                                     algorithm = cms.string("byDirectComparison"), # check that they
                                     srcTrack = cms.string("tracker"),             # have the same
                                     srcState = cms.string("atVertex"),            # tracker track
                                     matchedTrack = cms.string("tracker"),         # can't check ref
                                     matchedState = cms.string("atVertex"),        # because of the
                                     maxDeltaR        = cms.double(0.01),          # embedding.
                                     maxDeltaLocalPos = cms.double(0.01),
                                     maxDeltaPtRel    = cms.double(0.01),
                                     sortBy           = cms.string("deltaR"),
                                     )

process.tkToMediumMuonsNoRPC = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("promptTrackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("MediumMuonsNoRPC"), # to all global muons
                                     algorithm = cms.string("byDirectComparison"), # check that they
                                     srcTrack = cms.string("tracker"),             # have the same
                                     srcState = cms.string("atVertex"),            # tracker track
                                     matchedTrack = cms.string("tracker"),         # can't check ref
                                     matchedState = cms.string("atVertex"),        # because of the
                                     maxDeltaR        = cms.double(0.01),          # embedding.
                                     maxDeltaLocalPos = cms.double(0.01),
                                     maxDeltaPtRel    = cms.double(0.01),
                                     sortBy           = cms.string("deltaR"),
                                     )

process.tkToTightMuons = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("promptTrackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("TightMuons"), # to all global muons
                                     algorithm = cms.string("byDirectComparison"), # check that they
                                     srcTrack = cms.string("tracker"),             # have the same
                                     srcState = cms.string("atVertex"),            # tracker track
                                     matchedTrack = cms.string("tracker"),         # can't check ref
                                     matchedState = cms.string("atVertex"),        # because of the
                                     maxDeltaR        = cms.double(0.01),          # embedding.
                                     maxDeltaLocalPos = cms.double(0.01),
                                     maxDeltaPtRel    = cms.double(0.01),
                                     sortBy           = cms.string("deltaR"),
                                     )

process.tkToTightMuonsNoRPC = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("promptTrackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("TightMuonsNoRPC"), # to all global muons
                                     algorithm = cms.string("byDirectComparison"), # check that they
                                     srcTrack = cms.string("tracker"),             # have the same
                                     srcState = cms.string("atVertex"),            # tracker track
                                     matchedTrack = cms.string("tracker"),         # can't check ref
                                     matchedState = cms.string("atVertex"),        # because of the
                                     maxDeltaR        = cms.double(0.01),          # embedding.
                                     maxDeltaLocalPos = cms.double(0.01),
                                     maxDeltaPtRel    = cms.double(0.01),
                                     sortBy           = cms.string("deltaR"),
                                     )


process.passingMediumMuons = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToMediumMuons"),
                                       )

process.passingMediumMuonsNoRPC = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToMediumMuonsNoRPC"),
                                       )

process.passingTightMuons = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToTightMuons"),
                                       )

process.passingTightMuonsNoRPC = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToTightMuonsNoRPC"),
                                       )


## Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trackProbes@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200"),
)
#process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
#    decay = cms.string("tagMuons@+ probeMuons@-"), # charge coniugate states are implied
#    cut   = cms.string("40 < mass < 200"),
#)

## Match muons to MC
process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("tagMuons"), #"muons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

process.muMcMatchProbe = cms.EDProducer("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("trackProbes"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

## Make the tree
process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # pairs
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # variables to use
    variables = cms.PSet(
        ## methods of reco::Candidate
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ProbeCand = cms.InputTag("trackProbes"),
        PassingProbeMuons = cms.InputTag("probeMuons"),
        PassingMediumMuons = cms.InputTag("passingMediumMuons"),
        PassingTightMuons = cms.InputTag("passingTightMuons"),
        PassingMediumMuonsNoRPC = cms.InputTag("passingMediumMuonsNoRPC"),
        passingTightMuonsNoRPC = cms.InputTag("passingTightMuonsNoRPC"),
        ## two defined by simple string cuts
        #passingGlb = cms.string("isGlobalMuon"),
        #passingIso = cms.string("(isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),
    ),
    # mc-truth info
    isMC = cms.bool(MC_flag),
    motherPdgId = cms.vint32(22,23),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatchProbe"),
    allProbes     = cms.InputTag("trackProbes"),

)
##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
if MC_flag:
    process.tagAndProbe = cms.Path(
        process.primaryVertexFilter *
        process.HLTFilter *
        (process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
        process.promptMuons *
        process.promptMuonsNoRPC *
        process.PassingHLT *
        process.tightMuons *
        (process.tagMuons + process.probeMuons) *
        process.goodTracks *
        process.trackCands *
        process.promptTrackCands *
        process.trackProbes *
        process.MediumMuons *
        process.MediumMuonsNoRPC *
        process.TightMuons *
        process.TightMuonsNoRPC *
        process.tkToMediumMuons *
        process.tkToMediumMuonsNoRPC *
        process.tkToTightMuons *
        process.tkToTightMuonsNoRPC *
        process.passingMediumMuons *
        process.passingMediumMuonsNoRPC *
        process.passingTightMuons *
        process.passingTightMuonsNoRPC *
        (process.tpPairs + process.muMcMatch + process.muMcMatchProbe) *
        process.muonEffs
    )
else:
    process.tagAndProbe = cms.Path(
        process.primaryVertexFilter *
        process.HLTFilter *
        #process.runfilter *
        (process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
        process.promptMuons *
        process.promptMuonsNoRPC *
        process.PassingHLT *
        process.tightMuons *    
        (process.tagMuons + process.probeMuons) *
        process.goodTracks *
        process.trackCands *
        process.promptTrackCands *
        process.trackProbes *
        process.MediumMuons *
        process.MediumMuonsNoRPC *
        process.TightMuons *
        process.TightMuonsNoRPC *
        process.tkToMediumMuons *
        process.tkToMediumMuonsNoRPC *
        process.tkToTightMuons *
        process.tkToTightMuonsNoRPC *
        process.passingMediumMuons *
        process.passingMediumMuonsNoRPC *
        process.passingTightMuons *
        process.passingTightMuonsNoRPC *
        process.tpPairs *
        process.muonEffs
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TagProbeFitTree.root")
                                   )
