import FWCore.ParameterSet.Config as cms
process = cms.Process("RPCAnal")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "START42_V12::All"
process.GlobalTag.globaltag = "GR_R_42_V14::All"
process.prefer("GlobalTag")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/tmp/hkseo/Skim-ProbeF-Run2011A-PromptReco-v6_1.root',
    'file:/tmp/hkseo/Skim-ProbeF-Run2011A-PromptReco-v6_2.root',
    'file:/tmp/hkseo/Skim-ProbeF-Run2011A-PromptReco-v6_3.root'
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

### HLT filter
HLTPath = "HLT_IsoMu24_v*"
HLTProcessName = "HLT"
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZMuHLTFilter = copy.deepcopy(hltHighLevel)
process.ZMuHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","",HLTProcessName)
process.ZMuHLTFilter.throw = cms.bool(False)
process.ZMuHLTFilter.HLTPaths = [HLTPath]

process.HLTFilter = cms.Sequence(process.ZMuHLTFilter)


process.GlobalMuonsNoRPC = cms.EDFilter("MuonSelector",
                                        src = cms.InputTag("muonsNoRPC"),
                                        cut = cms.string("isGlobalMuon"
                                                         ),
                                        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('/tmp/hkseo/RPCRecHit_data_tree.root')
                                   )

process.demo = cms.EDAnalyzer('RecHitAnal',
                              Debug = cms.untracked.bool(False),
                              RPCRecHits = cms.InputTag("rpcRecHits"),
                              DT4DSegments = cms.InputTag("dt4DSegments"),
                              CSCSegments = cms.InputTag("cscSegments"),
                              muons = cms.InputTag("tightMuons"),
                              globalMuonsNoRPC = cms.InputTag("tightMuonsNoRPC"),
                              StandAloneMuons = cms.InputTag("standAloneMuonsNoRPC"),
                              generalTracks = cms.InputTag("generalTracks"),
                              SegmentsTrackAssociatorParameters = cms.PSet(
    segmentsDt = cms.untracked.InputTag("dt4DSegments"),
    SelectedSegments = cms.untracked.InputTag("SelectedSegments"),
    segmentsCSC = cms.untracked.InputTag("cscSegments")
    ),
                              )

# WZMu Skim
process.load("DPGAnalysis.Skims.WZMuSkim_cff")
process.Zskim = cms.Sequence(
    #process.WZMuHLTFilter *
    process.goodMuonsForZ *
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
    #process.Zskim *
    #process.primaryVertexFilter *
    #(process.muontrackingNoRPC+process.muonIdProducerSequenceNoRPC) *
    #process.GlobalMuonsNoRPC *
    process.demo
    )
