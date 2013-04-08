import FWCore.ParameterSet.Config as cms

kshortVertex = cms.EDFilter("VertexCandProducer",
    track = cms.PSet(
        src = cms.InputTag("generalTracks"),
        minPt = cms.double(3.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        signif = cms.double(0.5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        dxy = cms.double(0.0),
        signif = cms.double(5.0),
    ),
    pdgId = cms.uint32(310),
    leg1Id = cms.uint32(211),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(0.4),
    rawMassMax = cms.double(0.6),
    massMin = cms.double(0.43),
    massMax = cms.double(0.56),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(100),
)

lambdaVertex = kshortVertex.clone(
    pdgId = cms.uint32(3122),
    leg1Id = cms.uint32(2212),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(1.06),
    rawMassMax = cms.double(1.22),
    massMin = cms.double(1.08),
    massMax = cms.double(1.20),
)

phiVertex = kshortVertex.clone(
    pdgId = cms.uint32(333),
    leg1Id = cms.uint32(321),
    leg2Id = cms.uint32(321),
    rawMassMin = cms.double(0.98),
    rawMassMax = cms.double(1.06),
    massMin = cms.double(1.00),
    massMax = cms.double(1.04),
)
phiVertex.track.signif = -5
phiVertex.vertex.signif = -5
phiVertex.vertex.dxy = -999

jpsiVertex = kshortVertex.clone(
    pdgId = cms.uint32(443),
    leg1Id = cms.uint32(13),
    leg2Id = cms.uint32(13),
    rawMassMin = cms.double(2.75),
    rawMassMax = cms.double(3.45),
    massMin = cms.double(2.80),
    massMax = cms.double(3.40),
)
jpsiVertex.track.signif = -5
jpsiVertex.vertex.signif = -5
jpsiVertex.vertex.dxy = -999
