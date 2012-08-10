import FWCore.ParameterSet.Config as cms

KinematicTauAnalyzer = cms.EDAnalyzer("KinematicTauAnalyzer",
                                      discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality"),
                                      KinematicFitTauTag = cms.InputTag("KinematicTauProducer","KinematicFitTau"),
                                      gensrc = cms.InputTag('genParticles'),
                                      GenEventInfo  = cms.InputTag('generator'),
                                      TauMatchingDR = cms.double(0.2)
                                      )
