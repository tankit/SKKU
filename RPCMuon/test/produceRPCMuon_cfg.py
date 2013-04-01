import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("TEST")

#process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START61_V8::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/1661335F-3041-E211-9B96-00237DDBE0E2.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/2A12A045-2241-E211-8BF5-001D09F2915A.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/3C1D83B8-3641-E211-8C66-0025B32035A2.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/3E76F7E8-2741-E211-8249-003048D37666.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/3E863C44-2241-E211-9255-001D09F25041.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/4EB6D745-2241-E211-9738-001D09F24D8A.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/604E8D2C-2741-E211-B542-003048F11C28.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/6440884D-2941-E211-BBA9-0025901D6288.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/82B885ED-2241-E211-9877-001D09F252E9.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/8AAFC294-2141-E211-89E8-003048F1182E.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/90F6F479-2641-E211-99E5-001D09F29524.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/96BE2949-2241-E211-9993-001D09F23F2A.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/98EEEB5E-4A41-E211-A591-001D09F25460.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/A8CF653C-4D41-E211-811E-003048673374.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/AA4018D3-2C41-E211-8279-00215AEDFD98.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/AC6EF0B7-4941-E211-9EFB-003048D374F2.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/B27AC385-3241-E211-AD10-0019B9F4A1D7.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/E4E6B318-2041-E211-B351-001D09F29114.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/F2BA6B22-2C41-E211-9D7A-003048D2BED6.root',
        '/store/data/Run2012D/MinimumBias/RECO/PromptReco-v1/000/208/686/F60495B3-1E41-E211-BB7C-003048D3756A.root',
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

def lumiList( json ):
    import FWCore.PythonUtilities.LumiList as LumiList
    myLumis = LumiList.LumiList(filename = json ).getCMSSWString().split(',')
    return myLumis

def applyJSON( process, json ):

    # import PhysicsTools.PythonAnalysis.LumiList as LumiList
    # import FWCore.ParameterSet.Types as CfgTypes
    # myLumis = LumiList.LumiList(filename = json ).getCMSSWString().split(',')

    myLumis = lumiList( json )

    import FWCore.ParameterSet.Types as CfgTypes
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)

    # print process.source.lumisToProcess
jsonDir = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt"
jsonFile = "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt"
applyJSON(process, os.path.join(jsonDir, jsonFile))

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
#process.generalV0Candidates.lambdaMassCut = 0.10
#process.generalV0Candidates.kShortMassCut = 0.14

## User analysis
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root"),
)

process.load("SKKU.RPCMuon.VertexCandProducer_cfi")
process.load("SKKU.RPCMuon.fakeMuonAnalyzer_cfi")

process.muonRereco = cms.Sequence(
#    process.RawToDigi
    process.muonrecoComplete
  * process.muoncosmicreco * process.regionalCosmicTracksSeq
  * process.muoncosmichighlevelreco #* process.muonshighlevelreco
  * process.muons
  #* process.generalV0Candidates
)

process.pKs = cms.Path(
    process.primaryVertexFilter
  * process.kshortVertex
  * process.muonRereco
  * process.fakeKshort
)

process.pLambda = cms.Path(
    process.primaryVertexFilter
  * process.lambdaVertex
  * process.muonRereco
  * process.fakeLambda
)

process.pPhi = cms.Path(
    process.primaryVertexFilter
  * process.phiVertex
  * process.muonRereco
  * process.fakePhi
)

#process.outPath = cms.EndPath(process.out)

