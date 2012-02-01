import FWCore.ParameterSet.Config as cms
from RecoTauTag.KinematicTau.ThreeProngInputSelector_Step1_cfi import *
from RecoTauTag.KinematicTau.ThreeProngInputSelector_Step2_cfi import *
from RecoTauTag.KinematicTau.PrimaryVertexProducer_cfi import *

ThreeProngInputSelector = cms.Sequence(ThreeProngInputSelectorStep1*reducedPrimaryVertices*ThreeProngInputSelectorStep2)
ignoreThreeProngInputSelector = cms.Sequence(cms.ignore(ThreeProngInputSelectorStep1)
                                             *cms.ignore(reducedPrimaryVertices)
                                             *cms.ignore(ThreeProngInputSelectorStep2)
                                             )
