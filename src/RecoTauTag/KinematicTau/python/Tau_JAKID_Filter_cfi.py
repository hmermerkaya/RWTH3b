import FWCore.ParameterSet.Config as cms

TauJAKIDFilter = cms.EDFilter("Tau_JAKID_Filter",
                              jakid = cms.vint32(5),
                              nprongs = cms.vint32(3),
                              gensrc = cms.InputTag('genParticles'),
                              TauPtMin = cms.double(20.0),
                              TauEtaMax = cms.double(2.0)
                              )
