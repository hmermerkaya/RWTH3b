// -*- C++ -*-
//
// Package:    InputTrackSelector
// Class:      InputTrackSelector
// 
/**
 
 Description: creates collection of reco::TrackRefVector for each tau candidate
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Dec  15 19:21:54 CEST 2009
// $Id: InputTrackSelector.h,v 1.7 2010/07/12 13:07:35 perchall Exp $
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/KinematicFit/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/PFTauFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"


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
	std::string tauType_;
	edm::InputTag primVtx_;
	unsigned int minTracks_, minTau_, cnt_, cntFound_;
	bool filterTaus_;

};
