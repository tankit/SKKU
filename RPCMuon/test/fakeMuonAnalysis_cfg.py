import FWCore.ParameterSet.Config as cms
process = cms.Process("Ana")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.GlobalTag.globaltag = 'START53_V9::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1),)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/96CBE8F4-78F5-E111-99CB-003048FFD71A.root",
        "/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/E2C3D691-A4F5-E111-9624-0018F3D0968A.root",
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root"),
)

process.load("SKKU.RPCMuon.fakeMuonAnalyzer_cfi")

process.p = cms.Path(
    process.fakeKshort + process.fakeLambda
)
