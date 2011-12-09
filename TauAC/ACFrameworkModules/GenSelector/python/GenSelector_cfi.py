import FWCore.ParameterSet.Config as cms

GenSelector = cms.EDFilter("GenSelector",#creates reco::CandidateRefVector containing refs to selected genParticles
	candCollection = cms.InputTag("genParticles"),#genParticleCandidates #needs same collection as MCTruthDeltaRMatcherNew
	decayType = cms.untracked.string('VBFH'),#tauGun,Ztautau,...
    motherPdgID = cms.untracked.int32(0)#search for taus within a decay of a mother with given pdgID
)

printTree = cms.EDAnalyzer("ParticleListDrawer",
	maxEventsToPrint = cms.untracked.int32(0),
	printVertex = cms.untracked.bool(True),
	src = cms.InputTag("genParticles")
)

printSelectedTree = cms.EDAnalyzer("ParticleListDrawer",
	maxEventsToPrint = cms.untracked.int32(0),
	printVertex = cms.untracked.bool(True),
	src = cms.InputTag("GenSelector","genSignalDecay")
)

generatorSequence = cms.Sequence(printTree+cms.ignore(GenSelector))
#generatorSequence = cms.Sequence(cms.ignore(GenSelector)*printSelectedTree)#+ as printTree is not required
#generatorSequence = cms.Sequence(cms.ignore(GenSelector))
