import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *


KinematicTauProducer = cms.EDProducer("KinematicTauProducer",#creates reco::CandidateRefVector containing refs to selected jets
                                      #parameters for KinematicConstrainedVertexFitter
                                      fitParameters = cms.PSet( maxDelta = cms.double(.0001),#stopping condition
                                                                maxNbrOfIterations = cms.int32(100),	#number of iterations
                                                                maxReducedChiSq = cms.double(225.),
                                                                minChiSqImprovement = cms.double(50.)
                                                                ),
                                      KinematicTauCandTag = cms.InputTag("ThreeProngInputSelectorStep2","PreKinematicDecaysStep2") # the pre-selected taudecays 
                                      )
