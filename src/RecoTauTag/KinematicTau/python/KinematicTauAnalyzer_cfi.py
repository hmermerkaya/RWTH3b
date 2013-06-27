import FWCore.ParameterSet.Config as cms

KinematicTauAnalyzer = cms.EDAnalyzer("KinematicTauAnalyzer",
	kinematicTaus =	cms.InputTag("KinematicTauProducer"),
	discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality")
)
