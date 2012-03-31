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


ControlSample_Skim = cms.EDFilter('ControlSample_SkimCuts',
                                  muonsTag   = cms.InputTag("muons"),
                                  pfjetsTag  = cms.InputTag("ak5PFJets"),
                                  MuonPtCut  = cms.double(24.0),
                                  MuonEtaCut = cms.double(2.1),
                                  JetEtaCut  = cms.double(2.4),
                                  JetPtCut  = cms.double(10),
                                  dphiMuJet  = cms.double(2.5),
                                  )

    
