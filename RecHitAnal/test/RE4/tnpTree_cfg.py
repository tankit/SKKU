import FWCore.ParameterSet.Config as cms

GLOBAL_TAG = 'START61_V11::All'

HLTProcessName = "HLT"

process = cms.Process("TagProbe")

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Good vertex requirement
process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)
for i in range(1,100):
    #process.source.fileNames.append('/store/user/mskim/RE4/CMSSW_6_1_2-PRE_PO61_V1/20130628/reco_%d.root' % i)
    process.source.fileNames.append('/store/user/mskim/RE4/CMSSW_6_1_2-START61/20130629/reco_%d.root' % i)

### HLT filter
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = False
#process.hltHighLevel.HLTPaths = ["HLT_Mu24_v*", "HLT_IsoMu24_v*"]

##############
# TAG MUON
##############

process.promptMuons = cms.EDFilter("PromptMuonSelector",
    src = cms.InputTag("muons"),
    maxDxy = cms.untracked.double(0.2),
    maxDz = cms.untracked.double(0.2),
    #beamSpot = cms.InputTag("offlineBeamSpot"),
    vertex = cms.InputTag("goodPrimaryVertices"),
)

#for version in range(1,21):
#    process.PassingHLT.hltTags.append(cms.InputTag("HLT_IsoMu24_v%d::%s" % (version, HLTProcessName)))
#    process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu24_v%d" % version)
## Tags. In a real analysis we should require that the tag muon fires the trigger,
##       that's easy with PAT muons but not RECO/AOD ones, so we won't do it here
##       (the J/Psi example shows it)
#PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % ("HLT_Mu30_v7",);

process.tightMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("promptMuons"),
    cut = cms.string(
        "isGlobalMuon && isPFMuon"# && isolationR03().sumPt<3.0"
        "&& globalTrack().normalizedChi2() < 10.0"
        "&& globalTrack().hitPattern().numberOfValidMuonHits() > 0"
        "&& numberOfMatchedStations() > 1"
        #"&& abs(muonBestTrack().dxy(vertex.position())) < 0.2"
        #"&& abs(muonBestTrack().dz(vertex.position())) < 0.5"
        "&& innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        "&& track().hitPattern().trackerLayersWithMeasurement() > 5"
        "&& pt > 20 && abs(eta) < 2.4"
        #"&& track().hitPattern().numberOfValidPixelHits() > 0" #--added by Minsuk on Feb 17, 2012
        #"&& track().hitPattern().numberOfValidTrackerHits() > 10"
        #"&& innerTrack().numberOfValidHits()>10 && globalTrack().normalizedChi2()<10.0"
        #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
        #"&& numberOfMatchedStations>1" #--updated from numberOfMatches by Minsuk on May 12, 2012
        #"&& (isolationR03().sumPt+isolationR03().emEt+isolationR03().hadEt)<0.1*pt"
    ),
)

process.tagMuons = process.tightMuons.clone()

#process.tagMuons = cms.EDProducer("trgMatchedMuonProducer",
#    InputProducer = cms.InputTag("tightMuons"),
#    hltTags = cms.VInputTag(),
#    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
#    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
#)

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
    maxDxy = cms.untracked.double(0.2),
    maxDz = cms.untracked.double(0.2),
    #beamSpot = cms.InputTag("offlineBeamSpot"),
    vertex = cms.InputTag("goodPrimaryVertices"),
)

process.trackProbes = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("promptTrackCands"),
    cut = cms.string("abs(eta)<2.4"),
)

############
# MUON ID
############

process.LooseMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("promptMuons"),
    cut = cms.string(
        "(isGlobalMuon || isTrackerMuon) && isPFMuon"
        #"&& pt > 20 && abs(eta) < 2.4"
    ),
)

process.TightMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("LooseMuons"),
    cut = cms.string(
        "isPFMuon && isGlobalMuon"
        "&& globalTrack.normalizedChi2<10.0"
        "&& numberOfMatchedStations > 1"
        "&& globalTrack.hitPattern.numberOfValidMuonHits>0"
        "&& track.hitPattern.numberOfValidPixelHits > 0"
        "&& track.hitPattern.trackerLayersWithMeasurement > 5"
    ),
)

process.RPCMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("promptMuons"),
    cut = cms.string(
        "isRPCMuon"
        "&& numberOfMatchedRPCLayers >= 2"
        #"&& numberOfMatchedStations(RPCHitAndTrackArbitration) >= 2"
    ),
)

########################
# MATCH TRACK AND MUONS
########################

process.tkToLooseMuons = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("promptTrackCands"),
    matched = cms.InputTag("LooseMuons"),
    algorithm = cms.string("byDirectComparison"),
    srcTrack = cms.string("tracker"),
    srcState = cms.string("atVertex"),
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(0.01),
    maxDeltaLocalPos = cms.double(0.01),
    maxDeltaPtRel    = cms.double(0.01),
    sortBy           = cms.string("deltaR"),
)

process.tkToTightMuons = process.tkToLooseMuons.clone(
    matched = cms.InputTag("TightMuons"),
)

process.tkToRPCMuons = process.tkToLooseMuons.clone(
    matched = cms.InputTag("RPCMuons"),
)

process.passingLooseMuons = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("trackProbes"),
    match = cms.InputTag("tkToLooseMuons"),
)

process.passingTightMuons = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("trackProbes"),
    match = cms.InputTag("tkToTightMuons"),
)

process.passingRPCMuons = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("trackProbes"),
    match = cms.InputTag("tkToRPCMuons"),
)

## Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trackProbes@-"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
)

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
        abseta = cms.string("abs(eta)"),
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ProbeCand = cms.InputTag("trackProbes"),
        PassingLooseMuons = cms.InputTag("passingLooseMuons"),
        PassingTightMuons = cms.InputTag("passingTightMuons"),
        PassingRPCMuons = cms.InputTag("passingRPCMuons"),
        ## two defined by simple string cuts
        #passingGlb = cms.string("isGlobalMuon"),
        #passingIso = cms.string("(isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),
    ),
    # mc-truth info
    isMC = cms.bool(True),
    motherPdgId = cms.vint32(22,23),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatchProbe"),
    allProbes     = cms.InputTag("trackProbes"),

)

process.tagAndProbe = cms.Path(
    process.goodPrimaryVertices
#    process.hltHighLevel
  * process.promptMuons
  * process.tightMuons
  * process.tagMuons
  * process.goodTracks * process.trackCands * process.promptTrackCands * process.trackProbes
  * process.LooseMuons  * process.tkToLooseMuons  * process.passingLooseMuons
  * process.TightMuons  * process.tkToTightMuons  * process.passingTightMuons
  * process.RPCMuons    * process.tkToRPCMuons    * process.passingRPCMuons
  * (process.tpPairs + process.muMcMatch + process.muMcMatchProbe)
  * process.muonEffs
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("tnpTree.root")
)
