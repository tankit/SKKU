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
        #"/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/96CBE8F4-78F5-E111-99CB-003048FFD71A.root",
        #"/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/E2C3D691-A4F5-E111-9624-0018F3D0968A.root",
        "file:///tmp/jhgoh/96CBE8F4-78F5-E111-99CB-003048FFD71A.root",
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

process.p = cms.Path(
#    process.RawToDigi
    process.primaryVertexFilter
  * process.muonrecoComplete
  * process.muoncosmicreco * process.regionalCosmicTracksSeq
  * process.muoncosmichighlevelreco #* process.muonshighlevelreco
  * process.muons
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

