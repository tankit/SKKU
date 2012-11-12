import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

#process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V9::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/0201F024-2E1C-E211-A01D-001D09F24D4E.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/30CDBB47-2C1C-E211-B393-003048D37456.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/367B6D59-2D1C-E211-BDB8-5404A63886CC.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/5EDE6642-2E1C-E211-9078-001D09F25267.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/8A22545B-301C-E211-B060-003048D3733E.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/9AD9B459-2D1C-E211-BF0D-BCAEC518FF6E.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/B63124BD-2C1C-E211-AF8A-003048D2BC42.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/D2F1D925-2D1C-E211-A35F-003048678110.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/D4EDB704-2C1C-E211-8279-003048D2BBF0.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/E09F86AD-2D1C-E211-811D-0025901D5D78.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/EC865E73-2E1C-E211-AE10-001D09F25267.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/F668A9AF-2E1C-E211-9693-003048D373AE.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/F6B42A37-481C-E211-8C2C-003048D2BC4C.root',
        '/store/relval/CMSSW_6_1_0_pre4/RelValTTbar/GEN-SIM-RECO/PU_START61_V1-v3/0003/FC467103-2F1C-E211-A5B6-001D09F25267.root',
    ),
    inputCommands = cms.untracked.vstring(
        "drop *_*_*_RECO",

        "keep *_generalTracks_*_*",
        "keep *_tevMuons_*_*",
        "keep *_siPixel*_*_*", "keep *_siStrip*_*_*",
        "keep *_dt*_*_*", "keep *_csc*_*_*", "keep *_rpc*_*_*",
        "keep *_*Digi_*_*", "keep *_*Digis_*_*", 
        "keep *EcalRecHit*_*_*_*", "keep *H*RecHit*_h*reco_*_*",
        "keep *CaloTower*_*_*_*",
        "keep *_*CaloJets_*_*",
        "keep *_offlineBeamSpot_*_*",
        "keep *_offlinePrimaryVertices_*_*",
    ),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',),
)

from Configuration.EventContent.EventContent_cff import RECOSIMEventContent
process.out.outputCommands += RECOSIMEventContent.outputCommands

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

process.muidRPCMuLoose = cms.EDProducer("MuonSelectionTypeValueMapProducer",
    inputMuonCollection = cms.InputTag("muons1stStep"),
    selectionType = cms.string('RPCMuLoose'),
)

process.muonSelectionTypeSequence += process.muidRPCMuLoose

## Disable PFMuon for temporary solution not to crash in 53X
process.muons.FillPFIsolation = False
process.muons.FillPFMomentumAndAssociation = False

## Change Lambda and Kshort mass range
process.generalV0Candidates.lambdaMassCut = 0.05
process.generalV0Candidates.kShortMassCut = 0.07

process.p = cms.Path(
#    process.RawToDigi
    process.primaryVertexFilter
  * process.muonrecoComplete
  * process.muoncosmicreco * process.regionalCosmicTracksSeq
  * process.muoncosmichighlevelreco #* process.muonshighlevelreco
  * process.muons
  * process.generalV0Candidates
#    process.localreco * process.globalreco
#  + process.egammaHighLevelRecoPrePF + process.particleFlowReco
#  + process.regionalCosmicTracksSeq * process.muoncosmichighlevelreco * process.muonshighlevelreco
#  * process.particleFlowLinks
#  * process.jetHighLevelReco * process.tautagging
#  + process.metrecoPlusHCALNoise + process.btagging * process.recoPFMET + process.PFTau
#  * process.reducedRecHits
#  * process.reconstruction
)

process.outPath = cms.EndPath(process.out)

### User analysis

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('hist.root'),
#)

process.rpcMuAna = cms.EDAnalyzer("RPCMuonAnalyzer",
    muon = cms.untracked.InputTag("muons"),
    minPtTrk = cms.untracked.double(3),
    maxEtaTrk = cms.untracked.double(1.6),
)

