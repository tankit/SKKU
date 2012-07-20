import FWCore.ParameterSet.Config as cms

MC_flag = False
GLOBAL_TAG = 'GR_R_52_V8::All'
if MC_flag:
    GLOBAL_TAG = 'START52_V5::All'

process = cms.Process("TEST")

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = GLOBAL_TAG

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)   

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/m/mskim/public/rpc/temp/00AD4245-A5B5-E111-A1E8-001EC9D8B54A.root',
        'file:/afs/cern.ch/work/m/mskim/public/rpc/temp/0215F1C2-BBAE-E111-A01F-485B39800B69.root',
    ),
    inputCommands = cms.untracked.vstring(
    ),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',),
)

from Configuration.EventContent.EventContent_cff import RECOSIMEventContent
process.out.outputCommands += RECOSIMEventContent.outputCommands

process.muidRPCMuLoose = cms.EDProducer("MuonSelectionTypeValueMapProducer",
    inputMuonCollection = cms.InputTag("muons1stStep"),
    selectionType = cms.string('RPCMuLoose'),
)

process.muonSelectionTypeSequence += process.muidRPCMuLoose
process.p = cms.Path(
    process.RawToDigi
  * process.localreco * process.globalreco
  + process.egammaHighLevelRecoPrePF + process.particleFlowReco
  + process.regionalCosmicTracksSeq * process.muoncosmichighlevelreco * process.muonshighlevelreco
#  * process.particleFlowLinks
#  * process.jetHighLevelReco * process.tautagging
#  + process.metrecoPlusHCALNoise + process.btagging * process.recoPFMET + process.PFTau
#  * process.reducedRecHits
#  * process.reconstruction
)

process.outPath = cms.EndPath(process.out)

### User analysis

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist.root'),
)

process.rpcMuAna = cms.EDAnalyzer("RPCMuonAnalyzer",
    muon = cms.untracked.InputTag("muons"),
    minPtTrk = cms.untracked.double(3),
    maxEtaTrk = cms.untracked.double(1.6),
)

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

process.p += process.primaryVertexFilter*process.rpcMuAna
