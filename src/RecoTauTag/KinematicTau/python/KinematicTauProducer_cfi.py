import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.Geometry.GeometryIdeal_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
import os, subprocess


KinematicTauProducer = cms.EDProducer("KinematicTauProducer",#creates reco::CandidateRefVector containing refs to selected jets
                                      #parameters for KinematicConstrainedVertexFitter
                                      primVtx = cms.InputTag("offlinePrimaryVertices"),
                                      KinematicTauCandTag = cms.InputTag("ThreeProngInputSelectorStep1","PreKinematicDecaysStep1"),
                                      tauDaughterTracks = cms.InputTag("generalTracks"),
                                      beamSpot = cms.InputTag("offlineBeamSpot"),
                                      fitParameters = cms.PSet( maxDelta = cms.double(.01),#stopping condition
                                                                maxNbrOfIterations = cms.int32(100),	#number of iterations
                                                                maxReducedChiSq = cms.double(225.),
                                                                minChiSqImprovement = cms.double(50.)
                                                                ),
                                      gensrc = cms.InputTag('genParticles'),
                                      minVtxTracks = cms.untracked.int32(3),
                                      maxChi2ndf = cms.untracked.double(10.0),
                                      BDTweightFileMinus = cms.untracked.string(os.getenv("CMSSW_BASE")+"/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"), # currently not trained so all ambiguities are the same
                                      BDTweightFilePlus  = cms.untracked.string(os.getenv("CMSSW_BASE")+"/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"),
                                      BDTweightFileZero  = cms.untracked.string(os.getenv("CMSSW_BASE")+"/src/RecoTauTag/KinematicTau/QualityCutsTraining_BDT.weights.xml"),
                                      minTau = cms.untracked.uint32(1), #minimum taus to select (otherwise filter returns false)            
                                      )
