import FWCore.ParameterSet.Config as cms

TrigFilter = cms.EDFilter("Filter",
                          HLTResults = cms.InputTag( "TriggerResults" )


)
