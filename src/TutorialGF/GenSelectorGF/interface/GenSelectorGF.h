// -*- C++ -*-
//
// Package:    GenSelectorGF
// Class:      GenSelectorGF
// 
/**
 *
 * Description: <one line class summary>
 *
 * Implementation:
 * <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla
//         Created:  Thu Oct  8 11:12:54 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

class GenSelectorGF : public edm::EDFilter {
public:
	explicit GenSelectorGF(const edm::ParameterSet&);
	~GenSelectorGF();
	
private:
	virtual void beginJob() ;
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	bool checkGenEvt(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
	bool checkGenEvtZ3pr(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
	

	edm::Handle<reco::GenParticleCollection> genCandidate_;
	edm::InputTag candCollection_;
	std::string decay_;
	unsigned int cnt_, cntFound_;
	
};
