import FWCore.ParameterSet.Config as cms

GenSelector = cms.EDFilter("GenSelectorGF",#creates reco::CandidateRefVector containing refs to selected genParticles
	candCollection = cms.InputTag("genParticles"),#genParticleCandidates #needs same collection as MCTruthDeltaRMatcherNew
	decayType = cms.untracked.string('Z3pr')
)
