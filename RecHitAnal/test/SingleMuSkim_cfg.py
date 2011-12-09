import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Skim")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'START42_V12::All' #Summer11 MC
process.GlobalTag.globaltag = 'GR_R_42_V14::All' #2011A data

##from Configuration.PyReleaseValidation.autoCond import autoCond
##--This file is getting obsolete, please use Configuration.AlCa.autoCond instead
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

## Source
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#    pickRelValInputFiles( cmsswVersion  = os.environ['CMSSW_VERSION']
#                        , relVal        = 'RelValTTbar'
#                        , globalTag     = process.GlobalTag.globaltag.value().split(':')[0]
#                        , numberOfFiles = 1
#                        )
#    )
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:/korserv001/hkseo/data/Run2011A/SingleMu-ZMuSkim/Run2011A-ZMuSkim_55_1_0IG.root'
    'rfio:/castor/cern.ch/user/h/hkseo/data/Run2011A/ZMuSkim/Run2011A-ZMuSkim_55_1_0IG.root'
    )
)

# TagMu Skim
#process.load("SKKU.RecHitAnal.TagMuSkim_cff")
process.load("TagMuSkim_cff")

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )


process.p1 = cms.Path(process.primaryVertexFilter*process.TagMuSelSeq)
#process.p1 = cms.Path(process.diMuonSelSeq*process.demo)
#process.p2 = cms.Path(process.pfMetWMuNuSeq*process.demo)

# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('Skim.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p1') )
                               )
process.outpath = cms.EndPath(process.out)
