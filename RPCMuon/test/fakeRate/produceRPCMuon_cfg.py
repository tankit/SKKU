import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("TEST")

#process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'PRE_62_V8::All'
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_6_2_0/RelValProdTTbar/GEN-SIM-RECO/PRE_ST62_V8-v3/00000/0450C239-84EC-E211-A57D-003048F1DB62.root",
        "/store/relval/CMSSW_6_2_0/RelValProdTTbar/GEN-SIM-RECO/PRE_ST62_V8-v3/00000/402C3BE8-9EEC-E211-A668-C860001BD87E.root",
        "/store/relval/CMSSW_6_2_0/RelValProdTTbar/GEN-SIM-RECO/PRE_ST62_V8-v3/00000/76A79AFD-5AEC-E211-AD34-0025B3203804.root",
    ),
    inputCommands = cms.untracked.vstring(
        "drop *_*_*_RECO",

        "keep *_generalTracks_*_*",
        #"keep *_generalV0Candidates_*_*",
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

## Muon rereconstruction
#muidRPCMuLoose = cms.EDProducer("MuonSelectionTypeValueMapProducer",
#    inputMuonCollection = cms.InputTag("muons1stStep"),
#    selectionType = cms.string('RPCMuLoose'),
#)
#muonSelectionTypeSequence += muidRPCMuLoose

process.muonRereco = cms.Sequence(
#    process.RawToDigi
    process.muonrecoComplete
  * process.muoncosmicreco * process.regionalCosmicTracksSeq
  * process.muoncosmichighlevelreco #* process.muonshighlevelreco
  * process.muons
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

## Disable PFMuon for temporary solution not to crash in 53X
process.muons.FillPFIsolation = False
process.muons.FillPFMomentumAndAssociation = False

## User analysis
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root"),
)

process.load("SKKU.RPCMuon.VertexCandProducer_cfi")
process.load("SKKU.RPCMuon.fakeMuonAnalyzer_cfi")

process.pKs = cms.Path(
    process.noscraping + process.primaryVertexFilter
  + process.kshortVertex
  + process.muonRereco
  * process.fakeKshort
)

process.pLambda = cms.Path(
    process.noscraping + process.primaryVertexFilter
  + process.lambdaVertex
  + process.muonRereco
  #* process.fakeLambda2
  * process.fakeLambda
)

process.pPhi = cms.Path(
    process.noscraping + process.primaryVertexFilter
  + process.phiVertex
  + process.muonRereco
  * process.fakePhi
)

process.pJpsi = cms.Path(
    process.noscraping + process.primaryVertexFilter
  + process.jpsiVertex
  + process.muonRereco
  * process.fakeJpsi
)

#from SKKU.RPCMuon.applyJSON_cff import *
#jsonDir = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt"
#jsonFile = "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt"
##jsonFile = "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt"
#applyJSON(process, os.path.join(jsonDir, jsonFile))


#process.outPath = cms.EndPath(process.out)
if 'SECTION' in os.environ:
    section = int(os.environ['SECTION'])
    nFiles = 10

    jsonDir = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt"
    #jsonFile = "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt"
    jsonFile = "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt"
    applyJSON(process, os.path.join(jsonDir, jsonFile))

    files = [l.strip() for l in open("MinimumBias-Run2012D.txt").readlines()]
    begin = nFiles*section
    end = min(nFiles*(section+1), len(files)+1)
    process.source.fileNames = files[begin:end]

    process.TFileService.fileName = "unmerged/result_%03d.root" % section

    print process.source.fileNames[0]
    print process.source.fileNames[-1]

