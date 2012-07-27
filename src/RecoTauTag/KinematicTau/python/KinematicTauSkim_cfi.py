import FWCore.ParameterSet.Config as cms

KinematicTauSkim = cms.EDFilter("KinematicTauSkim",
                                kinematicTaus =	cms.InputTag("KinematicTauBasicProducer"),
                                #discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality")
                                discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFitQuality"),
                                KinematicFitTauTag = cms.InputTag("KinematicTauProducer","KinematicFitTauDecays"),
                                minTau = cms.untracked.uint32(1)#minimum taus to produce (otherwise filter returns false)	
                                )
