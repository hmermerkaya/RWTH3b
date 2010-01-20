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
// $Id: KinematicTauProducer.h,v 1.3 2010/01/19 09:27:22 perchall Exp $
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
#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"//own class of tauGroup to store fitted particles in event stream
#include "DataFormats/KinematicFit/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/PFTauFwd.h"


class KinematicTauProducer : public edm::EDFilter {
public:
	typedef std::vector<SelectedKinematicParticleCollection> KinematicCollection;
	explicit KinematicTauProducer(const edm::ParameterSet&);
	~KinematicTauProducer();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	bool select(KinematicCollection & refitParticles, const reco::Vertex & primVtx, InputTauCollection & PFTauRef, reco::PFCandidateCollection & PFDaughters);
//	int addRefittedParticles(const int &ambiguityCnt, RefCountedKinematicTree tree, KinematicConstrainedVertexFitter* kcvFitter, KinematicCollection &refitParticles, reco::PFTauRef &tauRef, const reco::Vertex &primVtx);
//	void correctReferences(KinematicCollection & selected, edm::OrphanHandle<reco::PFCandidateCollection> & orphanPFCands);
	bool combineHiggs(reco::Vertex & primVtx, SelectedKinematicParticleCollection & selectedHiggs);
	
	const edm::ParameterSet& iConfig_;
	edm::Event * iEvent_;
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
	
	edm::InputTag primVtx_, usedTauCandidatesTag_, inputCollectionTag_;
	unsigned int cnt_, cntFound_;
	
	std::vector<RefCountedKinematicTree> *tauM_, *tauP_;
	
};
