import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

KinematicTauProducer = cms.EDFilter("KinematicTauProducer",#creates reco::CandidateRefVector containing refs to selected jets
	#parameters for KinematicConstrainedVertexFitter
	fitParameters = cms.PSet(#parameters for KinematicConstrainedVertexFitter
		maxDelta = cms.double(.001),#stopping condition
		maxNbrOfIterations = cms.int32(20),	#number of iterations
		maxReducedChiSq = cms.double(225.),
		minChiSqImprovement = cms.double(50.)
	),
	primVtx = cms.InputTag("ThreeProngInputSelector","primVtx"),#selected offlinePrimaryVerticesFromCTFTrack, use the reduced vertex from ThreeProngInputSelector here
	selectedTauCandidates = cms.InputTag("ThreeProngInputSelector","InputTauRefs"),
	inputTracks = cms.InputTag("ThreeProngInputSelector","InputTracks"),#selected tracks from PFTaus (daughters of selectedTauCandidates)
)
