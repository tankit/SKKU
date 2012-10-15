import FWCore.ParameterSet.Config as cms
import os

section = int(os.environ['SECTION'])
begin = int(os.environ['BEGIN'])
end = int(os.environ['END'])
sample = os.environ['SAMPLE']

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.GlobalTag.globaltag = "START52_V12::All"

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
files = []
for line in open('samples/%s.txt' % sample).readlines():
    line = line.strip()
    files.append(line)
process.source.fileNames.extend(files[begin:end])

#print process.source.fileNames[0]
#print process.source.fileNames[-1]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result/result_%s_%03d.root" % (sample, section)),
)

process.genParticleCount = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(13, -13, 11, -11),
    minNumber = cms.uint32(1), 
    maxNumber = cms.uint32(1),
)

process.genParticleTauVeto = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(15, -15),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25),
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", 
    src = cms.InputTag('offlinePrimaryVertices'),
    filterParams =  cms.PSet(
        minNdof = cms.double(4.),
        maxZ    = cms.double(24.), 
        maxRho  = cms.double(2.)
    ),
    filter = cms.bool(True),
)

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
HLTPaths = {
    "MuHad_5E33":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v*"],
    "MuHad_7E33":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v*"],

    "EleHad_5E33":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*"],
    "EleHad_7E33":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet45_35_25_v*"],

    "Mu_51X":["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4",],
    "Mu_52X_GTV5":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2",],
    "Mu_52X_GTV9":["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v2",],

    "Ele_51X":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v4",],
    "Ele_52X_GTV5":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8",],
    "Ele_52X_GTV9":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v3",],
}
#process.hltHighLevel.HLTPaths = HLTPaths["MuHad_5E33"] + HLTPaths["MuHad_7E33"]
process.hltHighLevel.HLTPaths = HLTPaths["Mu_52X_GTV9"]

process.event = cms.EDAnalyzer("EventTupleProducerMuon",
    doMCMatch = cms.bool(True),
    gen = cms.InputTag("genParticles"),
    jet = cms.InputTag("loosePatJetsPF"),
    met = cms.InputTag("patMETsPF"),
    lepton = cms.InputTag("goodPatMuonsPF"),
    #lepton = cms.InputTag("goodElectronsPF"),
    leptonCut = cms.string(
        "abs(eta) < 2.1 && pt > 30 && dB < 0.2"
        " && isPFMuon && isGlobalMuon && normChi2 < 10"
        " && track.hitPattern.trackerLayersWithMeasurement > 5"
        " && globalTrack.hitPattern.numberOfValidMuonHits > 0"
        " && innerTrack.hitPattern.numberOfValidPixelHits > 0"
        " && numberOfMatchedStations > 1"
        " && (chargedHadronIso+neutralHadronIso+photonIso) < 0.2*pt"
    ),
    jetCut = cms.string(
      " abs(eta) < 2.5 && pt > 35 && numberOfDaughters() > 1"
      " && neutralHadronEnergyFraction() < 0.99 && neutralEmEnergyFraction() < 0.99"
      " && (abs(eta) >= 2.4 || chargedEmEnergyFraction() < 0.99)"
      " && (abs(eta) >= 2.4 || chargedHadronEnergyFraction() > 0)"
      " && (abs(eta) >= 2.4 || chargedMultiplicity() > 0)"
    ),
    bTagType = cms.string("combinedSecondaryVertexBJetTags"),
)

process.p = cms.Path(
  #  process.genParticleCount + process.genParticleTauVeto
  #+ process.printDecay
#  +
 process.goodOfflinePrimaryVertices
#  + process.hltHighLevel
  * process.event
)
