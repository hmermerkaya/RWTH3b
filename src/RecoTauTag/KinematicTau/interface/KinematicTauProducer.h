// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      KinematicTauProducer
// 
/**
 * This framework module copies an existing PFTauCollection and modifies the tau's parameters according to a kinematical refit of its decay.
 * A particular decay mode has to be assumed. New discriminators regarding the fits quality are provided.
 * Part of the KinematicTau package.
 *
 * @author Lars Perchalla, Philip Sauerland
 * @date 2009
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Dec  16 11:12:54 CEST 2009
// $Id: KinematicTauProducer.h,v 1.18 2010/08/13 12:29:00 perchall Exp $
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
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"


class KinematicTauProducer : public edm::EDFilter {
public:
	explicit KinematicTauProducer(const edm::ParameterSet&);
	~KinematicTauProducer();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	/**
	 Execute the kinematic fit and in case of success modify the taus parameters
	 */
	bool select(reco::PFTauCollection & selected, std::map<int, std::vector<bool> > & discrimValues, const reco::Vertex & primaryVtx);
	/**
	 combine a discriminator of important quality cuts of a refitted tau decay
	 */
	bool dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr, const int & fitStatus, const reco::PFTauRef & tauRef, const reco::Vertex & primaryVtx);
	/**
	 Fill the new PFTau dicriminators from fit result
	 */
	void discriminate(const edm::OrphanHandle<reco::PFTauCollection> & collection, const std::map<int, std::vector<bool> > & dicrimValues);
	/**
	 Calulate the maximal allowed GJ angle from a1 parameters
	 */	
	double thetaGJMax(double ma1, double pa1, double Mtau = 1.777);
	
	const edm::ParameterSet fitParameters_;
	edm::Event * iEvent_;
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
	
	edm::InputTag primVtx_, selectedTauCandidatesTag_, inputCollectionTag_;
	
};
