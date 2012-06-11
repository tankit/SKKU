import FWCore.ParameterSet.Config as cms
##                      _              _
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##
################################################
MC_flag = False
#MC_flag = True
GLOBAL_TAG = 'GR_R_42_V21A::All'
if MC_flag:
    GLOBAL_TAG = 'START42_V12::All'

#HLTPath1 = "HLT_IsoMu24_v1"
#HLTPath2 = "HLT_IsoMu24_v8"
#HLTPath3 = "HLT_IsoMu24_v9"
HLTProcessName = "HLT"

RECOProcess = "RECO"
JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.6 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99 && pt>15.0"

process = cms.Process("TagProbe")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
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
        'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/0201AF20-BFDB-E011-A808-0019B9F72F97.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/346EC6A2-8BDB-E011-8437-003048F11C5C.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/409BF93C-8ADB-E011-BEE9-003048F117EA.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/42512260-6DDB-E011-944D-001D09F231C9.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/48BC67A0-A3DB-E011-B31F-BCAEC532972E.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/5CBF4E7E-90DB-E011-B235-BCAEC5329709.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/8EC621FE-D6DB-E011-A210-001D09F2512C.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/A044F736-70DB-E011-82A8-BCAEC518FF91.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/A6BE5182-90DB-E011-A54E-BCAEC518FF3C.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/A8DD1F32-7EDB-E011-A951-003048F118C6.root',
        #'file:/korserv001/ehkwon/SingleMu_2011B_RECO_PromptReco-v1_175834/BA65A25D-6DDB-E011-9DB9-003048F118AC.root'
        
        #'rfio:/castor/cern.ch/user/e/ehkwon/Skim/MC/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/SingleMuSkim_9_1_iKA.root',
        #'rfio:/castor/cern.ch/user/e/ehkwon/Skim/MC/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/SingleMuSkim_2150_1_Tmw.root'
        
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_1.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_2.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_3.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_4.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_5.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_6.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_7.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_8.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_9.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_10.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_11.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_12.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_13.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_14.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_15.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_16.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_17.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_18.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_19.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_20.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_21.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_22.root',
        #'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_23.root'
        
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_0.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_1.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_2.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_3.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_4.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_5.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_6.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_7.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_8.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_9.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_10.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_11.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_12.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_13.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_14.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_15.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_16.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_17.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_18.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_19.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_20.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_21.root',
        #'/../user/e/ehkwon/Skim/2012Jan20_MC/ZMuSkim_22.root'
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

for version in range(1,11):
    process.PassingHLT.hltTags.append(cms.InputTag("HLT_IsoMu24_v%d::%s" % (version, HLTProcessName)))
    process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu24_v%d" % version)
    
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




