import FWCore.ParameterSet.Config as cms

TauJAKIDFilter = cms.EDFilter("Tau_JAKID_Filter",
                              jakid = cms.vint32(5,18,14),
                              gensrc = cms.InputTag('genParticles'),
                              TauPtMin = cms.double(15.0),
                              TauEtaMax = cms.double(2.0)
                              )
