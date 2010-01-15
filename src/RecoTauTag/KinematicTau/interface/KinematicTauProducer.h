// -*- C++ -*-
//
// Package:    KinematicTauProducer
// Class:      KinematicTauProducer
// 
/**\class KinematicTauProducer KinematicTauProducer.cc KinHiggs/KinematicTauProducer/src/KinematicTauProducer.cc
 
 Description: <one line class summary>
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla
//         Created:  Thu Dec  16 11:12:54 CEST 2009
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

class KinematicTauProducer : public edm::EDFilter {
public:
	typedef std::vector<SelectedKinematicParticleCollection> KinematicCollection;
	explicit KinematicTauProducer(const edm::ParameterSet&);
	~KinematicTauProducer();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	bool checkPrimVtx(reco::Vertex &primVtx);	


};
