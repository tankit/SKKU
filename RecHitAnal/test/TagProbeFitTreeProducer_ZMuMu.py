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
##process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')
#process.GlobalTag.globaltag = cms.string('GR_R_42_V21A::All')
#process.GlobalTag.globaltag = cms.string('START42_V12::All')
process.GlobalTag.globaltag = GLOBAL_TAG
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
#        'file:/korserv001/ehkwon/MC/DYToMuMu_M60_RECO_Eta16_oldCLS_1.root'

        'file:/korserv001/ehkwon/Run2011A-PromptReco-v6/SingleMuSkim_1_1_tVc.root',
#        'file:/korserv001/ehkwon/Run2011A-PromptReco-v6/SingleMuSkim_100_1_r0U.root',
#        'file:/korserv001/ehkwon/Run2011A-PromptReco-v6/SingleMuSkim_101_1_smF.root',
#        'file:/korserv001/ehkwon/Run2011A-PromptReco-v6/SingleMuSkim_102_1_B8w.root',
#        'file:/korserv001/ehkwon/Run2011A-PromptReco-v6/SingleMuSkim_103_1_WID.root'

#        'file:/korserv001/mskim/rpc/SingleMu_2011B_RECO_PromptReco-v1_175834_0201AF20-BFDB-E011-A808-0019B9F72F97.root'

#'/../user/e/ehkwon/Skim/2012Jan08/Run2011A-PromptReco-v6/SingleMuSkim_100_1_r0U.root',
#'/../user/e/ehkwon/Skim/2012Jan08/Run2011A-PromptReco-v6/SingleMuSkim_101_1_smF.root',
#'/../user/e/ehkwon/Skim/2012Jan08/Run2011A-PromptReco-v6/SingleMuSkim_102_1_B8w.root',
#'/../user/e/ehkwon/Skim/2012Jan08/Run2011A-PromptReco-v6/SingleMuSkim_103_1_WID.root'
    )
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
process.load("RecoMuon.MuonIdentification.muons_cfi")
process.muonsNoRPC = process.muons.clone()
process.muonsNoRPC.inputCollectionLabels = cms.VInputTag(
    cms.InputTag("generalTracks"),
    cms.InputTag("globalMuonsNoRPC"),
    cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
    )
process.muonsNoRPC.fillGlobalTrackQuality = False

process.muonIdProducerSequenceNoRPC = cms.Sequence(
    process.muonsNoRPC
    )
###############################################            


##############
# TAG MUON
##############

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
    process.PassingHLT.hltTags.append(cms.InputTag("HLT_IsoMu24_v%d::%s" % (version, HLTProcessName)))
    process.ZMuHLTFilter.HLTPaths.append("HLT_IsoMu24_v%d" % version)

## Tags. In a real analysis we should require that the tag muon fires the trigger, 
##       that's easy with PAT muons but not RECO/AOD ones, so we won't do it here
##       (the J/Psi example shows it)
#PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % ("HLT_Mu30_v7",);

process.tightMuons = cms.EDFilter("MuonSelector",
                                src = cms.InputTag("muons"),
                                cut = cms.string("isGlobalMuon && isTrackerMuon && isolationR03().sumPt<3.0"
                                                 "&& abs(innerTrack().dxy)<1.0 && pt > 20 && abs(eta) < 2.4"
                                                 "&& track().hitPattern().numberOfValidPixelHits() > 0" #--added by Minsuk on Feb 17, 2012
                                                 "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                 "&& innerTrack().numberOfValidHits()>10 && globalTrack().normalizedChi2()<10.0"
                                                 "&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                 "&& numberOfMatches>1"
                                                 "&& (isolationR03().sumPt+isolationR03().emEt+isolationR03().hadEt)<0.1*pt"
                                                 ), 
                                )

process.tagMuons = process.PassingHLT.clone()
process.tagMuons.InputProducer = cms.InputTag("tightMuons")


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

process.trackProbes = cms.EDFilter("CandViewRefSelector",
                                   src = cms.InputTag("trackCands"),
                                   cut = cms.string("abs(eta)<1.8"),
                                   )


############
# MUON ID
############

process.LooseMuons = cms.EDFilter("MuonSelector",
                                  src = cms.InputTag("muons"),
                                  cut = cms.string("isGlobalMuon"
                                                   #"&& numberOfMatches > 1"
                                                   "&& abs(innerTrack().dxy)<2.0 && pt > 0 && abs(eta) < 2.4"
                                                   "&& globalTrack().normalizedChi2()<10.0"
                                                   #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                   "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                   "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                   "&& innerTrack().numberOfValidHits()>10"
                                                   ), 
                                  )


process.LooseMuonsNoRPC = cms.EDFilter("MuonSelector",
                                       src = cms.InputTag("muonsNoRPC"),
                                       cut = cms.string("isGlobalMuon"
                                                        #"&& numberOfMatches > 1"
                                                        "&& abs(innerTrack().dxy)<2.0 && pt > 0 && abs(eta) < 2.4"
                                                        "&& globalTrack().normalizedChi2()<10.0"
                                                        #"&& globalTrack().hitPattern().numberOfValidMuonHits()>0"
                                                        "&& track().hitPattern().numberOfValidPixelHits() > 0"
                                                        "&& track().hitPattern().numberOfValidTrackerHits() > 10"
                                                        "&& innerTrack().numberOfValidHits()>10"
                                                        ), 
                                       )


process.MediumMuons = cms.EDFilter("MuonSelector",
                                   src = cms.InputTag("LooseMuons"),
                                   cut = cms.string("numberOfMatches > 1"), 
                                   )


process.MediumMuonsNoRPC = cms.EDFilter("MuonSelector",
                                        src = cms.InputTag("LooseMuonsNoRPC"),
                                        cut = cms.string("numberOfMatches > 1"), 
                                        )



########################
# MATCH TRACK AND MUONS
########################

process.tkToLooseMuons = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("trackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("LooseMuons"), # to all global muons
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

process.tkToLooseMuonsNoRPC = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("trackCands"), # all tracks are available for matching
                                     matched = cms.InputTag("LooseMuonsNoRPC"), # to all global muons
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

process.tkToMediumMuons = cms.EDProducer("MatcherUsingTracks",
                                     src     = cms.InputTag("trackCands"), # all tracks are available for matching
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
                                     src     = cms.InputTag("trackCands"), # all tracks are available for matching
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


process.passingLooseMuons = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToLooseMuons"),
                                       )

process.passingLooseMuonsNoRPC = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToLooseMuonsNoRPC"),
                                       )

process.passingMediumMuons = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToMediumMuons"),
                                       )

process.passingMediumMuonsNoRPC = cms.EDProducer("MatchedCandidateSelector",
                                       src   = cms.InputTag("trackProbes"),
                                       match = cms.InputTag("tkToMediumMuonsNoRPC"),
                                       )


## Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trackProbes@-"), # charge coniugate states are implied
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
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ProbeCand = cms.InputTag("trackProbes"),
        PassingLooseMuons = cms.InputTag("passingLooseMuons"),
        PassingMediumMuons = cms.InputTag("passingMediumMuons"),
        PassingLooseMuonsNoRPC = cms.InputTag("passingLooseMuonsNoRPC"),
        passingMediumMuonsNoRPC = cms.InputTag("passingMediumMuonsNoRPC"),
        ## two defined by simple string cuts
        #passingGlb = cms.string("isGlobalMuon"),
        #passingIso = cms.string("(isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),
    ),
    # mc-truth info
    isMC = cms.bool(False), #--False in case of track as probe although it is MC (otherwise, an error!)
    motherPdgId = cms.vint32(22,23),
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
if MC_flag:
    process.tagAndProbe = cms.Path(
        process.HLTFilter *
        (process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
        process.PassingHLT *
        process.tightMuons *
        process.tagMuons *
        process.goodTracks *
        process.trackCands *
        process.trackProbes *
        process.LooseMuons *
        process.LooseMuonsNoRPC *
        process.MediumMuons *
        process.MediumMuonsNoRPC *
        process.tkToLooseMuons *
        process.tkToLooseMuonsNoRPC *
        process.tkToMediumMuons *
        process.tkToMediumMuonsNoRPC *
        process.passingLooseMuons *
        process.passingLooseMuonsNoRPC *
        process.passingMediumMuons *
        process.passingMediumMuonsNoRPC *
        (process.tpPairs + process.muMcMatch) *
        process.muonEffs
    )
else:
    process.tagAndProbe = cms.Path(
        process.HLTFilter *
        #process.runfilter *
        (process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
        process.PassingHLT *
        process.tightMuons *    
        process.tagMuons *
        process.goodTracks *
        process.trackCands *
        process.trackProbes *
        process.LooseMuons *
        process.LooseMuonsNoRPC *
        process.MediumMuons *
        process.MediumMuonsNoRPC *
        process.tkToLooseMuons *
        process.tkToLooseMuonsNoRPC *
        process.tkToMediumMuons *
        process.tkToMediumMuonsNoRPC *
        process.passingLooseMuons *
        process.passingLooseMuonsNoRPC *
        process.passingMediumMuons *
        process.passingMediumMuonsNoRPC *
        process.tpPairs *
        process.muonEffs
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TagProbeFitTree.root")
                                   )
