import FWCore.ParameterSet.Config as cms

PrimVtxSelector = cms.EDFilter("PrimVtxSelector",
	primVtx = cms.InputTag("offlinePrimaryVertices"), #Primary vertex reconstructed using the tracks taken from the generalTracks collection
#	primVtx = cms.InputTag("offlinePrimaryVerticesWithBS"),	#Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position
	minTracks = cms.untracked.int32(3),
	maxChi2ndf = cms.untracked.double(10.0)
)
