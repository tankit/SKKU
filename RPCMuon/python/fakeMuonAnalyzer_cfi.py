import FWCore.ParameterSet.Config as cms

fakeKshort = cms.EDAnalyzer("FakeMuonAnalyzer",
    muon = cms.InputTag("muons"),
    vertexCand = cms.InputTag("generalV0Candidates", "Kshort"),
    muonIds = cms.PSet(),
    vertexCut = cms.string(""),
    match = cms.string("matchByTrackRef"),
    #maxDR = cms.double(0.02),
    #maxDPt = cms.double(0.1),
)

fakeKshort.muonIds.all = cms.PSet(
    cut = cms.string(""),
    idSelection = cms.string("All"),
    arbitration = cms.string("NoArbitration"),
)

fakeKshort.muonIds.globalMuon = cms.PSet(
    cut = cms.string("isGlobalMuon"),
    idSelection = cms.string("AllGlobalMuons"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.globalMuonPT = cms.PSet(
    cut = cms.string("isGlobalMuon"),
    idSelection = cms.string("GlobalMuonPromptTight"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.RPCMuLoose = cms.PSet(
    cut = cms.string("isRPCMuon"),
    idSelection = cms.string("RPCMuLoose"),
    arbitration = cms.string("RPCHitAndTrackArbitration"),
)

fakeLambda = fakeKshort.clone()
fakeLambda.vertexCand = cms.InputTag("generalV0Candidates", "Lambda")

