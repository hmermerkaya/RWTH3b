import FWCore.ParameterSet.Config as cms

KinematicTauSkim = cms.EDAnalyzer("KinematicTauSkim",
	kinematicTaus =	cms.InputTag("KinematicTauProducer"),
	#discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality")
	discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFitQuality")
)
