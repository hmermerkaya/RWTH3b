import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *


KinematicTauProducer = cms.EDProducer("KinematicTauProducer",#creates reco::CandidateRefVector containing refs to selected jets
                                      #parameters for KinematicConstrainedVertexFitter
                                      primVtx = cms.InputTag("offlinePrimaryVertices"),
                                      KinematicTauCandTag = cms.InputTag("ThreeProngInputSelectorStep1","PreKinematicDecaysStep1"),
                                      fitParameters = cms.PSet( maxDelta = cms.double(.0001),#stopping condition
                                                                maxNbrOfIterations = cms.int32(100),	#number of iterations
                                                                maxReducedChiSq = cms.double(225.),
                                                                minChiSqImprovement = cms.double(50.)
                                                                ),
                                      gensrc = cms.InputTag('genParticles'),
                                      minVtxTracks = cms.untracked.int32(3),
                                      maxChi2ndf = cms.untracked.double(10.0),
                                      BDTweightFileMinus = cms.untracked.string("$CMSSW_BASE/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"),
                                      BDTweightFilePlus  = cms.untracked.string("$CMSSW_BASE/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"),
                                      BDTweightFileZero  = cms.untracked.string("$CMSSW_BASE/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"),
                                      minTau = cms.untracked.uint32(1), #minimum taus to select (otherwise filter returns false)            
                                      )
