import FWCore.ParameterSet.Config as cms

#pfTau = 'shrinkingCone'		#Signal cone: dR = 5/Et, Iso Cone = 0.5, 0.07 < dR < 0.15; where Et is the tau tranvserse energy (still up-to-date?)
pfTau = 'hps'
#pfTau = 'fixedCone'			#Signal cone: dR = 0.07, Iso Cone = 0.5
#pfTau = 'fixedConeHighEff'		#Signal cone: dR = 0.15, Iso Cone = 0.5


InputTrackSelector = cms.EDFilter("InputTrackSelector",#creates PFTauRefVector and collection of vector<reco::TrackRefVector> for each tau cand
	tauType = cms.untracked.string(pfTau),	#default is shrinkingCone
	minTracks = cms.uint32(3),	#only tau candidates with more/equal than minTracks are selected
	minTau = cms.untracked.uint32(1)	#minimum taus to select (otherwise filter returns false)
	
	#reimplement optional application of common PFDiscriminators. hint: cvs revert ;-)
)
