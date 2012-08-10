import FWCore.ParameterSet.Config as cms
from RecoTauTag.KinematicTau.InputTrackSelector_cfi import *
from RecoTauTag.KinematicTau.PrimaryVertexProducer_cfi import *
from RecoTauTag.KinematicTau.ThreeProngInputSelector_Step2_cfi import *
from RecoTauTag.KinematicTau.kinematictau_cfi import *
from RecoTauTag.KinematicTau.KinematicTauSkim_cfi import *
from RecoTauTag.KinematicTau.KinematicTauAnalyzer_cfi import *

#define sequences for Kinematic Fit
KinematicFitSequence         = cms.Sequence(InputTrackSelector*reducedPrimaryVertices*ThreeProngInputSelectorStep2*KinematicTauProducer)
KinematicFitSequencewithSkim = cms.Sequence(InputTrackSelector*reducedPrimaryVertices*ThreeProngInputSelectorStep2*KinematicTauProducer*KinematicTauSkim)
KinematicFitSequencewithDQM  = cms.Sequence(InputTrackSelector*reducedPrimaryVertices*ThreeProngInputSelectorStep2*KinematicTauProducer*KinematicTauAnalyzer)

