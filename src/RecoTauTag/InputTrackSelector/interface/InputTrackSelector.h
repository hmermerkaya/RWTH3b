// -*- C++ -*-
//
// Package:    InputTrackSelector
// Class:      InputTrackSelector
// 
/**\class InputTrackSelector InputTrackSelector.cc HiggsKinTau/InputTrackSelector/src/InputTrackSelector.cc
 
 Description: creates collection of vector<reco::CandidateRef> for each tau cand
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla
//         Created:  Thu Dec  15 19:21:54 CEST 2009
// $Id$
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

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
//#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"//for: get<reco::TrackRef>()
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/KinematicFit/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/PFTauFwd.h"


class InputTrackSelector : public edm::EDFilter {
public:
	explicit InputTrackSelector(const edm::ParameterSet&);
	~InputTrackSelector();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	bool select(InputTrackCollection & refitParticles, InputTauCollection & PFTauRef);
	bool filterInput(reco::PFTauRef &tau);
	reco::TrackRefVector getPFTauDaughters(reco::PFTauRef &PFTau);
	
	edm::Event * iEvent_;

	edm::InputTag inputCollectionTag_, primVtx_;
	int verbosity_;
	unsigned int cnt_, cntFound_;

};
