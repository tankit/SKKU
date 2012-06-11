import FWCore.ParameterSet.Config as cms
##                      _              _
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##
################################################
#MC_flag = False
MC_flag = True
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
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.GlobalTag.globaltag = GLOBAL_TAG

######### EXAMPLE CFG 
###  A simple test of runnning T&P on Zmumu to determine muon isolation and identification efficiencies
###  More a showcase of the tool than an actual physics example

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.hltTrigReport = cms.EDAnalyzer( 'HLTrigReport',
#    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
#)
#process.HLTAnalyzerEndpath = cms.EndPath( process.hltTrigReport )
#process.MessageLogger.categories.append("HLTrigReport")

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/0E187509-0D7A-E111-8FA3-001A928116C2.root',
        '/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/1011EE9E-2B7A-E111-9349-0018F3D0970C.root',
        '/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/5CAA0235-0F7A-E111-BA3E-0018F3D09690.root',
        '/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/A29B9025-0E7A-E111-97E7-001A928116DE.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
#process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")

### HLT filter
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZMuHLTFilter = copy.deepcopy(hltHighLevel)
process.ZMuHLTFilter.throw = cms.bool(False)
#process.ZMuHLTFilter.HLTPaths = [HLTPath1,HLTPath2,HLTPath3]

process.HLTFilter = cms.Sequence(process.ZMuHLTFilter)

# Trigger  ##################
process.PassingHLT = cms.EDProducer("trgMatchedMuonProducer",
                                    InputProducer = cms.InputTag("muons"),
                                    hltTags = cms.VInputTag(
                                    #cms.InputTag(HLTPath1,"", HLTProcessName),
                                    #cms.InputTag(HLTPath2,"", HLTProcessName),
                                    #cms.InputTag(HLTPath3,"", HLTProcessName)
                                    ),
                                    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
                                    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
                                    )

for version in range(1,21):
    process.PassingHLT.hltTags.append(cms.InputTag("HLT_IsoMu24_eta2p1_v%d::%s" % (version, HLTProcessName)))
    process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu24_eta2p1_v%d" % version)
    #process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu20_eta2p1_v%d" % version)
    #process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu30_eta2p1_v%d" % version)
    #process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu34_eta2p1_v%d" % version)
    #process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu40_eta2p1_v%d" % version)
    ##process.ZMuHLTFilter.HLTPaths.append("HLT_Mu24_eta2p1_v%d" % version)
    
## Tags. In a real analysis we should require that the tag muon fires the trigger, 
##       that's easy with PAT muons but not RECO/AOD ones, so we won't do it here
##       (the J/Psi example shows it)
process.glbMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isGlobalMuon && pt > 20 && abs(eta) < 2"), 
)

#if MC_flag:
#    process.tagMuons = process.glbMuons.clone()
#else:
#    process.tagMuons = process.PassingHLT.clone()
#    process.tagMuons.InputProducer = cms.InputTag("glbMuons")

process.tagMuons = process.PassingHLT.clone()
process.tagMuons.InputProducer = cms.InputTag("glbMuons")

## Probes. Now we just use Tracker Muons as probes
process.probeMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isTrackerMuon && pt > 10"), 
)

## Here we show how to define passing probes with a selector
## although for this case a string cut in the TagProbeFitTreeProducer would be enough
process.probesPassingCal = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string(process.probeMuons.cut.value() + " && caloCompatibility > 0.6"),
)

## Here we show how to use a module to compute an external variable
process.drToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("muons"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag("ak5CaloJets"),
    objectSelection = cms.InputTag("et > 20 && abs(eta) < 3 && n60 > 3 && (.05 < emEnergyFraction < .95)"),
)


## Here we show how to use a module to compute an external variable
#producer of dR < 0.4 muon-cleaned jets
process.cleanJets = cms.EDProducer("JetViewCleaner",
    srcObject = cms.InputTag(JET_COLL, "", "RECO"),
    srcObjectSelection = cms.string(JET_CUTS),
    srcObjectsToRemove = cms.VInputTag( cms.InputTag("muons", "", RECOProcess)),
    deltaRMin = cms.double(0.4)
    )


#produce dR(muon, nearest IDed uncorrected jet passing cuts on corrected eta and pT)
process.muonDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    #probes = cms.InputTag("probeMuons"),
    probes = cms.InputTag("muons"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS)
)

## Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200"),
)

## Match muons to MC
process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("muons"),
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
        pt  = cms.string("pt"),
        phi = cms.string("phi"),
        ## a method of the reco::Muon object (thanks to the 3.4.X StringParser)
        nsegm_old = cms.string("numberOfMatches"),
        nsegm = cms.string("numberOfMatchedStations"), 
        ## this one is an external variable
        #drj = cms.InputTag("drToNearestJet"),
        drj = cms.InputTag("muonDRToNearestJet"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ## one defined by an external collection of passing probes
        #passingCal = cms.InputTag("probesPassingCal"), 
        ## two defined by simple string cuts
        passingGlb = cms.string("isGlobalMuon && pt > 20 && abs(eta) < 1.8"),
        passingIso = cms.string("(isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),
        passingGlbIso = cms.string("isGlobalMuon && pt > 20 && abs(eta) < 1.8 && (isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),
    ),
    # mc-truth info
    isMC = cms.bool(MC_flag),
    motherPdgId = cms.vint32(22,23), #gamma, Z0
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("probeMuons"),
)
##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
## 'A*B' means 'B needs output of A'; # 'A+B' means 'if you want you can re-arrange the order'
if MC_flag:
    process.tagAndProbe = cms.Path(
        process.primaryVertexFilter *
        process.HLTFilter *
        process.PassingHLT *
        process.glbMuons *
        (process.tagMuons + process.probeMuons) *
        (process.cleanJets + process.muonDRToNearestJet +
         process.tpPairs + process.muMcMatch) *
        process.muonEffs
    )
else:
    process.tagAndProbe = cms.Path(
        process.primaryVertexFilter *
        process.HLTFilter *
        process.PassingHLT *
        process.glbMuons *
        (process.tagMuons + process.probeMuons) *
        (process.cleanJets + process.muonDRToNearestJet +
         process.tpPairs) *
        process.muonEffs
    )



process.TFileService = cms.Service("TFileService", fileName = cms.string("testTagProbeFitTreeProducer_ZMuMu.root"))




