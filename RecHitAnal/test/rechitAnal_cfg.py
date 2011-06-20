import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo2")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"
process.prefer("GlobalTag")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/../user/h/hkseo/RPC/ZMuSkim/zmumuJet_384_10_1_Ibq.root'
    #'file:/tmp/hkseo/FC960582-297F-E011-886A-001A92971B5E.root'
    'file:/tmp/hkseo/ZMuSkim-Summer11_9_1_3DX.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('/tmp/hkseo/RPCRecHit_tree.root')
                                   )

process.demo = cms.EDAnalyzer('RecHitAnal',
                              Debug = cms.untracked.bool(False),
                              RPCRecHits = cms.InputTag("rpcRecHits"),
                              DT4DSegments = cms.InputTag("dt4DSegments"),
                              CSCSegments = cms.InputTag("cscSegments"),
                              muons = cms.InputTag("muons"),
                              StandAloneMuons = cms.InputTag("standAloneMuons"),
                              generalTracks = cms.InputTag("generalTracks")
                              )


# WZMu Skim
#process.load("DPGAnalysis.Skims.WZMuSkim_cff")
#process.Zskim = cms.Sequence(
#    process.WZMuHLTFilter *
#    process.goodMuonsForZ *
#    process.dimuons *
#    process.dimuonsFilter
#    )

# Good vertex requirement
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                           minimumNDOF = cms.uint32(4) ,
#                                           maxAbsZ = cms.double(24),
#                                           maxd0 = cms.double(2)
#                                           )

process.p1 = cms.Path(process.demo)
#process.p1 = cms.Path(process.Zskim*process.demo)
#process.p1 = cms.Path(process.pfMetWMuNuSeq*process.primaryVertexFilter*process.demo)
#process.p2 = cms.Path(process.tcMetWMuNuSeq*process.primaryVertexFilter*process.demo)
