import FWCore.ParameterSet.Config as cms
process = cms.Process("Ana")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V9::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1),)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_5_3_2-START53_V6/RelValProdQCD_Pt_3000_3500/AODSIM/v2/0000/F01CD266-BDB9-E111-83EC-003048FFCBFC.root",
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root"),
)

process.allMuon = cms.EDAnalyzer("RPCFakeAnalyzer",
    muon = cms.InputTag("muons"),
    kshort = cms.InputTag("generalV0Candidates", "Kshort"),
    idCuts = cms.vstring(
        "isGlobalMuon",
        "isTrackerMuon",
    #    "dB < 0.2",
    )
)

process.p = cms.Path(
    process.allMuon
)
