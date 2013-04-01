import FWCore.ParameterSet.Config as cms

kshortVertex = cms.EDProducer("VertexCandProducer",
    track = cms.PSet(
        src = cms.InputTag("generalTracks"),
        minPt = cms.double(3.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        ipSignif = cms.double(2.),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        dxy = cms.double(0.0),
        vtxSignif = cms.double(15.0),
    ),
    pdgId = cms.uint32(310),
    leg1Id = cms.uint32(211),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(0.4),
    rawMassMax = cms.double(0.6),
    massMin = cms.double(0.43),
    massMax = cms.double(0.56),
)

lambdaVertex = kshortVertex.clone(
    pdgId = cms.uint32(2114),
    leg1Id = cms.uint32(2212),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(1.116-0.1),
    rawMassMax = cms.double(1.116+0.1),
    massMin = cms.double(1.116-0.06),
    massMax = cms.double(1.116+0.06),
)

phiVertex = kshortVertex.clone(
    pdgId = cms.uint32(3),
    leg1Id = cms.uint32(321),
    leg2Id = cms.uint32(321),
    rawMassMin = cms.double(1.020-0.1),
    rawMassMax = cms.double(1.020+0.1),
    massMin = cms.double(1.020-0.06),
    massMax = cms.double(1.020+0.06),
)
