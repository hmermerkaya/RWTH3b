import FWCore.ParameterSet.Config as cms

PreselectionCuts = cms.EDFilter('SkimmingCuts',
                    hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                    MuonPtCut = cms.double(24.0),
                    MuonIsGlobal = cms.bool(True),
                    NMuons      = cms.double(1.0),
                    MuonEtaCut = cms.double(2.1),
                    PFTauPtCut =cms.double(15.0),
                    PFTauEtaCut=cms.double(2.1)

)
