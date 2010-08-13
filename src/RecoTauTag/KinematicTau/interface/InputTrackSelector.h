// -*- C++ -*-
//
// Package:    InputTrackSelector
// Class:      InputTrackSelector
// 
/**
 * The InputTrackSelector selects tracks from within a PFTaus signal cone and stores combinations of these according to the required number of tracks defined by the user.
 *
 * This framework module returns a boolean whether the minimum number of taus was found with enough daughters.
 * For convenience a reference to the taus is also stored in case of success.
 *
 * @author Lars Perchalla, Philip Sauerland in 2009
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Dec  15 19:21:54 CEST 2009
// $Id: InputTrackSelector.h,v 1.8 2010/08/10 15:19:07 sauerlan Exp $
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
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"


class InputTrackSelector : public edm::EDFilter {
public:
	explicit InputTrackSelector(const edm::ParameterSet&);
	~InputTrackSelector();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	bool select(std::vector<reco::TrackRefVector> & refitParticles, reco::PFTauRefVector & PFTauRef);
	reco::TrackRefVector getPFTauDaughters(reco::PFTauRef &PFTau);
	
	edm::Event * iEvent_;
	std::string tauType_;
	edm::InputTag primVtx_;
	unsigned int minTracks_, minTau_, cnt_, cntFound_;

};
