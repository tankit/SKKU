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
        "/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/96CBE8F4-78F5-E111-99CB-003048FFD71A.root",
        "/store/relval/CMSSW_5_3_4_cand1-START53_V10/RelValProdQCD_Pt_3000_3500/GEN-SIM-RECO/v1/0000/E2C3D691-A4F5-E111-9624-0018F3D0968A.root",
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root"),
)

process.allMuon = cms.EDAnalyzer("FakeMuonAnalyzer",
    muon = cms.InputTag("muons"),
    vertexCand = cms.InputTag("generalV0Candidates", "Kshort"),
    muonTypes = cms.PSet(
        all = cms.string(""),
        TM = cms.string("isTrackerMuon"),
        GLB = cms.string("isGlobalMuon"),
        STA = cms.string("isStandAloneMuon"),
    ),
    vertexCut = cms.string(""),
    match = cms.string("matchByTrackRef"),
    #maxDR = cms.double(0.02),
    #maxDPt = cms.double(0.1),
)

process.p = cms.Path(
    process.allMuon
)
