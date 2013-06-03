// -*- C++ -*-
//
// Package:    KinematicTauProducer
// Class:      KinematicTauProducer
// 
/**
 
 Description: test application of KinematicTauCreator
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Dec  16 11:12:54 CEST 2009
// $Id: KinematicTauProducer.h,v 1.6 2010/01/26 13:15:07 perchall Exp $
//
//


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"//own class of tauGroup to store fitted particles in event stream
#include "DataFormats/KinematicFit/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/PFTauFwd.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


class KinematicTauProducer : public edm::EDFilter {
public:
	explicit KinematicTauProducer(const edm::ParameterSet&);
	~KinematicTauProducer();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	bool select(SelectedKinematicDecayCollection & refitDecays, InputTauCollection & PFTauRefCollection, reco::RecoChargedCandidateCollection & daughterCollection, const reco::Vertex & primaryVtx);
	void saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection);
	int saveKinParticles(KinematicTauCreator *kinTauCrtr, SelectedKinematicDecayCollection &refitDecays, const reco::Vertex &primVtx);
	void correctReferences(SelectedKinematicDecayCollection & selected, edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands);
	
	const edm::ParameterSet fitParameters_;
	edm::Event * iEvent_;
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
	
	edm::InputTag primVtx_, selectedTauCandidatesTag_, inputCollectionTag_;
	unsigned int minKinTau_, cnt_, cntFound_;
	
	static std::string intToString(int f){
		char s[32]; sprintf(s, "%d", f);
		return std::string(s);
	}	
};
