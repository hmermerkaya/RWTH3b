import FWCore.ParameterSet.Config as cms

PreselectionCuts = cms.EDFilter('SkimmingCuts',
                    hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                    MuonPtCut = cms.double(17.0),
                    MuonIsGlobal = cms.bool(True),
                    MuonEtaCut = cms.double(2.1),
                    PFTauPtCut =cms.double(17.0),
                    PFTauEtaCut=cms.double(2.1)

)
