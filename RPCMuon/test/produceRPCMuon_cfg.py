import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START52_V5::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)   

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/02B15C95-077A-E111-A28F-001A92971B9A.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/16A420C8-097A-E111-8EED-0026189438DD.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/2241659B-077A-E111-B8D0-002618943896.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/3448AD0D-077A-E111-B6F4-003048678FDE.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/58B2B90D-087A-E111-819D-00304867BFC6.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/78BDC1FE-067A-E111-A8BD-003048678BE8.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/AA292F02-077A-E111-9985-003048FFD71E.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/C2A09C49-2C7A-E111-9028-002618943943.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/C6B15699-087A-E111-A135-003048FFD760.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/E8446578-067A-E111-B897-002618943972.root',
        '/store/relval/CMSSW_5_2_3/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V5-v1/0043/FCAF0197-087A-E111-9854-0026189438D9.root',
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

process.muidRPCMuMedium = cms.EDProducer("MuonSelectionTypeValueMapProducer",
    inputMuonCollection = cms.InputTag("muons1stStep"),
    selectionType = cms.string('RPCMuMedium'),
)

process.muonSelectionTypeSequence += process.muidRPCMuMedium
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

#process.outPath = cms.EndPath(process.out)

### User analysis

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist.root'),
)

process.rpcMuAna = cms.EDAnalyzer("RPCMuonAnalyzer",
    muon = cms.untracked.InputTag("muons"),
)

process.p += process.rpcMuAna
