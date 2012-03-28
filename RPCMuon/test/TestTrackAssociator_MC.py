import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MC_43_V4::All' #Dima's test
process.GlobalTag.globaltag = 'START42_V12::All' #Summer11 MC
#process.GlobalTag.globaltag = 'GR_R_42_V14::All' #2011A data
#process.GlobalTag.globaltag = 'GR_R_42_V21A:All' #2011B data

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
# process.load("MagneticField.Engine.volumeBasedMagneticField_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("Geometry.DTGeometry.dtGeometry_cfi")

process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

process.load("Geometry.CSCGeometry.cscGeometry_cfi")

process.load("Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

#process.load("DQMServices.Components.MEtoEDMConverter_cfi")
#process.load("DQMServices.Core.DQM_cfg")

# add TrackDetectorAssociator lookup maps to the EventSetup
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff") 
# from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *  
from TrackingTools.TrackAssociator.default_cfi import * 

process.demo = cms.EDAnalyzer('TestTrackAssociator',
    TrackAssociatorParameterBlock,maxNEvents = cms.untracked.int32(-1),
                              outputAscii = cms.string("dummy.ascii"),
                              minPtTrk = cms.untracked.int32(20)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)  
)

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    # categories = cms.untracked.vstring('TrackAssociator','TrackAssociatorVerbose'),
    # categories = cms.untracked.vstring('TrackAssociator'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        # threshold = cms.untracked.string('DEBUG'),
	noTimeStamps = cms.untracked.bool(True),
	noLineBreaks = cms.untracked.bool(True),
	DEBUG = cms.untracked.PSet(
           limit = cms.untracked.int32(0)
	),
	# TrackAssociator = cms.untracked.PSet(
	#   limit = cms.untracked.int32(-1)
	#),
	# TrackAssociatorVerbose = cms.untracked.PSet(
	#   limit = cms.untracked.int32(-1)
	#),
    ),
    debugModules = cms.untracked.vstring("demo")
)

import FWCore.Framework.test.cmsExceptionsFatalOption_cff
options = cms.untracked.PSet(
    Rethrow = FWCore.Framework.test.cmsExceptionsFatalOption_cff.Rethrow
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'/store/relval/CMSSW_4_3_0_pre6/RelValSingleMuPt10/GEN-SIM-RECO/MC_43_V3-v1/0086/BA16DCB9-2A8C-E011-AD40-0030486791BA.root'
        #'file:~/666945BD-0E7E-E011-95F4-00237DA28240.root'
        #'file:/korserv001/mskim/rpc/Run2011A_SingleMu_WMu-May10ReReco-v1/666945BD-0E7E-E011-95F4-00237DA28240.root'

        #--/castor/cern.ch/user/h/hkseo/mc/Summer11/ZMuSkim

        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_100_1_PwS.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_101_1_614.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_103_1_gEP.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_104_1_Jws.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_105_1_KAj.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_106_1_ifO.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_108_1_dPu.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_109_1_ATC.root',
        #'file:/korserv001/hkseo/mc/Summer11/DYToMuMu-ZMuSkim/ZMuSkim-Summer11_10_1_VCj.root',

    )
)

for line in open('DYToMuMu-ZMuSkim-Summer11.txt').readlines():
    line = line.strip("'\", \n")
    if '.root' not in line: continue

    process.source.fileNames.append(line)

#process.DQMStore = cms.Service("DQMStore")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('output.root')
                                   )

process.p = cms.Path(process.demo)


