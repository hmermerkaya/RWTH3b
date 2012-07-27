import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

ThreeProngInputSelectorStep2 = cms.EDFilter("ThreeProngInputSelector_Step2",
                                            threeProngs = cms.InputTag("ThreeProngInputSelectorStep1", "ThreeProngCombinations"),
                                            primVtx = cms.InputTag("reducedPrimaryVertices"),
                                            selectedTauCandidates = cms.InputTag("InputTrackSelector", "InputTauRefs"),
                                            KinematicTauCandTag = cms.InputTag("InputTrackSelector","PreKinematicDecaysStep1"),
                                            minVtxTracks = cms.untracked.int32(3),
                                            maxChi2ndf = cms.untracked.double(10.0),
                                            minTau = cms.untracked.uint32(1), #minimum taus to select (otherwise filter returns false)
                                            )
