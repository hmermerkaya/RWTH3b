import FWCore.ParameterSet.Config as cms

TrigFilter = cms.EDFilter("TriggerFilter",
                          HLTResults = cms.InputTag( "TriggerResults" ),
                          doTauplusXTrigger = cms.untracked.bool(True),
                          doMuonTrigger = cms.untracked.bool(False),
                          doElectronTrigger = cms.untracked.bool(False),
                          doTauplusMETTrigger = cms.untracked.bool(False)
                          )

TrigFilterInfo = cms.EDProducer("TriggerFilterInfoProducer")
