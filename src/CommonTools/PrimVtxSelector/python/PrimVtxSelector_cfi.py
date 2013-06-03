import FWCore.ParameterSet.Config as cms

PrimVtxSelector = cms.EDFilter("PrimVtxSelector",
	primVtx = cms.InputTag("offlinePrimaryVertices"),#offlinePrimaryVerticesFromCTFTrack
	minTracks = cms.untracked.int32(3),
	maxChi2ndf = cms.untracked.double(10.0),
	verbosity = cms.untracked.int32(1)
)
