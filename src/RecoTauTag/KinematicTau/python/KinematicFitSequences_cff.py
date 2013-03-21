import FWCore.ParameterSet.Config as cms

from RecoTauTag.KinematicTau.ThreeProngInputSelector_Step1_cfi import *
from RecoTauTag.KinematicTau.KinematicTauProducer_cfi import *
from RecoTauTag.KinematicTau.KinematicTauSkim_cfi import *
from RecoTauTag.KinematicTau.KinematicTauAnalyzer_cfi import *
from RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi import *

#sequences for Kinematic Fit 
KinematicFitSequence         = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer)
KinematicFitSequencewithSkim = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer*KinematicTauSkim)
KinematicFitSequencewithDQM  = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer*KinematicTauAnalyzer)
