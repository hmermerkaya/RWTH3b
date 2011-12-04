import FWCore.ParameterSet.Config as cms

KinematicTauSkim = cms.EDFilter("KinematicTauSkim",
	kinematicTaus =	cms.InputTag("KinematicTauProducer"),
	#discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality")
	discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFitQuality"),
	minTau = cms.untracked.uint32(1)#minimum taus to produce (otherwise filter returns false)	
)