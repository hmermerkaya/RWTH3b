import FWCore.ParameterSet.Config as cms

ThreeProngInputSelectorStep1 = cms.EDFilter("ThreeProngInputSelector_Step1", 
	tauCandidates = cms.InputTag("InputTrackSelector", "InputTracks"),
    tracks = cms.InputTag("generalTracks"),
    minTau = cms.untracked.uint32(1), #minimum taus to select (otherwise filter returns false)
)
