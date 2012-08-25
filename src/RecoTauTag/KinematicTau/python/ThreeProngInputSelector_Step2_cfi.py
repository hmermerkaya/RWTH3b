import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

ThreeProngInputSelectorStep2 = cms.EDProducer("ThreeProngInputSelector_Step2",
                                              primVtx = cms.InputTag("offlinePrimaryVertices"),
                                              KinematicTauCandTag = cms.InputTag("ThreeProngInputSelectorStep1","PreKinematicDecaysStep1"),
                                              minVtxTracks = cms.untracked.int32(3),
                                              maxChi2ndf = cms.untracked.double(10.0),
                                              minTau = cms.untracked.uint32(1), #minimum taus to select (otherwise filter returns false)
                                            )
