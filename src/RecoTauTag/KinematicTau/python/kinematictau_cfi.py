import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

KinematicTauProducer = cms.EDFilter("KinematicTauProducer",#creates reco::CandidateRefVector containing refs to selected jets
	#parameters for KinematicConstrainedVertexFitter
	fitParameters = cms.PSet(#parameters for KinematicConstrainedVertexFitter
		maxDistance = cms.double(.001),#stopping condition
		maxNbrOfIterations = cms.int32(20),	#number of iterations
		maxOfInitialValue = cms.double(9999.)
	),
	primVtx = cms.InputTag("PrimVtxSelector","primVtx"),#selected offlinePrimaryVerticesFromCTFTrack
	selectedTauCandidates = cms.InputTag("InputTrackSelector","InputTauRefs"),
	inputTracks = cms.InputTag("InputTrackSelector","InputTracks"),#selected tracks from PFTaus (daughters of selectedTauCandidates)
)
