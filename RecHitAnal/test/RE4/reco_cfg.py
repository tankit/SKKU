# Auto generated configuration file
# using: 
# Revision: 1.13 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --datatier GEN-SIM-RECO --conditions PRE_PO61_V1::All -s RAW2DIGI,L1Reco,RECO -n -1 --eventcontent FEVTSIM --geometry ExtendedPostLS1 --filein file:step1.root --fileout reco.root --no_exec --mc
import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process('RECO')

simFile = os.environ['SIMFILE']
recoFile = os.environ['RECOFILE']

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryExtendedPostLS1Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:%s' % simFile)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.13 $'),
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTSIMEventContent.outputCommands,
    fileName = cms.untracked.string(recoFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition
process.FEVTSIMoutput.outputCommands.extend([
    'keep CSCDetIdCSCRPCDigiMuonDigiCollection_muonCSCDigis_MuonCSCRPCDigi_*',
    'keep RPCDetIdRPCDigiMuonDigiCollection_muonRPCDigis_*_*',
])


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_PO61_V1::All', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'START61_V11::All', '')
from CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi import *
process.myCond = EcalTrivialConditionRetriever.clone()
process.es_prefer_gt = cms.ESPrefer("PoolDBESSource","GlobalTag")

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTSIMoutput_step = cms.EndPath(process.FEVTSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTSIMoutput_step)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()
# Turn off some flags for CSCRecHitD that are turned ON in default config
process.csc2DRecHits.readBadChannels = cms.bool(False)
process.csc2DRecHits.CSCUseGasGainCorrection = cms.bool(False)
