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
// $Id: KinematicTauProducer.h,v 1.10 2010/06/09 15:49:22 perchall Exp $
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
#include "DataFormats/KinematicFit/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
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
	
	bool select(reco::PFTauCollection & selected, std::map<int, std::vector<bool> > & discrimValues, const reco::Vertex & primaryVtx);
	bool dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr);
	void discriminate(const edm::OrphanHandle<reco::PFTauCollection> & collection, const std::map<int, std::vector<bool> > & dicrimValues);
	
	const edm::ParameterSet fitParameters_;
	edm::Event * iEvent_;
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
	
	edm::InputTag primVtx_, selectedTauCandidatesTag_, inputCollectionTag_;
	
};
