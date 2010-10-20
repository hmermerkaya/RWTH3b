import FWCore.ParameterSet.Config as cms

TutorialAnlzr = cms.EDAnalyzer("TutorialAnlzr",
	generator =	cms.InputTag("GenSelector","genSignalDecay")
)
