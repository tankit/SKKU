import FWCore.ParameterSet.Config as cms
process = cms.Process("RPCAnal")

##                      _              _
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##
################################################
MC_flag = False
#MC_flag = True
GLOBAL_TAG = 'GR_R_53_V3::All'
if MC_flag:
    GLOBAL_TAG = 'START53_V9::All'
    
HLTPath = "HLT_IsoMu24*"
HLTProcessName = "HLT"
    

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.GlobalTag.globaltag = GLOBAL_TAG

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/tmp/hkseo/skim_singleMu_run193207_7_1_CPm.root'
    )
)

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

# muons1stStep
process.load("RecoMuon.MuonIdentification.muons1stStep_cfi")
process.muonsNoRPC1stStep = process.muons1stStep.clone()
process.muonsNoRPC1stStep.inputCollectionTypes = cms.vstring('inner tracks','links','outer tracks')
process.muonsNoRPC1stStep.inputCollectionLabels = cms.VInputTag(
    cms.InputTag("generalTracks"),
    cms.InputTag("globalMuonsNoRPC"),
    cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
    )
process.muonsNoRPC1stStep.fillGlobalTrackQuality = cms.bool(False)

# muons
process.load("RecoMuon.MuonIdentification.muons2muons_cfi")
process.muons.FillPFIsolation = False
process.muons.FillSelectorMaps = False
process.muons.FillCosmicsIdMap =  False
process.muons.FillTimingInfo = False
process.muons.FillDetectorBasedIsolation = False
process.muons.FillShoweringInfo = False

process.muonsNoRPC = process.muons.clone()
process.muonsNoRPC.InputMuons = cms.InputTag("muonsNoRPC1stStep")
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
###############################################

### HLT filter
#import copy
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#process.ZMuHLTFilter = copy.deepcopy(hltHighLevel)
#process.ZMuHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","",HLTProcessName)
#process.ZMuHLTFilter.throw = cms.bool(False)
#process.ZMuHLTFilter.HLTPaths = [HLTPath]

#process.HLTFilter = cms.Sequence(process.ZMuHLTFilter)


process.GlobalMuonsNoRPC = cms.EDFilter("MuonSelector",
                                        src = cms.InputTag("muonsNoRPC"),
                                        cut = cms.string("isGlobalMuon"
                                                         ),
                                        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('RpcHits_Tree.root')
                                   )

process.RpcHits = cms.EDAnalyzer('RecHitAnal',
                              Debug = cms.untracked.bool(False),
                              RPCRecHits = cms.InputTag("rpcRecHits"),
                              DT4DSegments = cms.InputTag("dt4DSegments"),
                              CSCSegments = cms.InputTag("cscSegments"),
                              muons = cms.InputTag("tightMuonsForZ"),
                              globalMuonsNoRPC = cms.InputTag("muonsNoRPC"),
                              StandAloneMuons = cms.InputTag("standAloneMuonsNoRPC"),
                              generalTracks = cms.InputTag("generalTracks"),
                              SegmentsTrackAssociatorParameters = cms.PSet(
    segmentsDt = cms.untracked.InputTag("dt4DSegments"),
    SelectedSegments = cms.untracked.InputTag("SelectedSegments"),
    segmentsCSC = cms.untracked.InputTag("cscSegments")
    ),
                              )

# WZMu Skim
process.load("DPGAnalysis.Skims.ZMuSkim_cff")
process.Zskim = cms.Sequence(
    process.ZMuHLTFilter *
    process.looseMuonsForZ *
    process.tightMuonsForZ *
    process.dimuons *
    process.dimuonsFilter
    )

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

process.p1 = cms.Path(
    #process.HLTFilter *
    process.primaryVertexFilter *
    process.Zskim *
    (process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
    process.GlobalMuonsNoRPC *
    process.RpcHits
    )
