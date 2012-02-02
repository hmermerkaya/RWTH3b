// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      ThreeProngInputSelector_Step1
// 
/**
 * This framework module creates the appropriate input for the KinematicTauProducer.
 * In this case the tau decay tau->3pi+nu is assumed.
 * It creates a collection of reco::TrackRefVector of size 3 for each tau candidate and refits the primary vertex without tracks assigned already to the tau decay.
 * Part of the KinematicTau package.
 *
 * @author Lars Perchalla, Philip Sauerland
 * @date 2012
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Feb  18 13:25:12 CEST 2010
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TLorentzVector.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef std::vector<std::vector<std::vector<reco::TrackRef> > > vVVTrackRef;
typedef std::vector<std::vector<reco::TrackRef> > vVTrackRef;

class ThreeProngInputSelector_Step1 : public edm::EDFilter {
public:
	explicit ThreeProngInputSelector_Step1(const edm::ParameterSet &);
	~ThreeProngInputSelector_Step1();
	
private:
	virtual void beginJob();
	virtual bool filter(edm::Event &, const edm::EventSetup &);
	virtual void endJob();
	bool sumCharge(const std::vector<reco::TrackRef> & input);
    template <typename T> std::vector<std::vector<T> > permuteCombinations(const std::vector<T> & vect);
    vVTrackRef choose3Prongs(std::vector<reco::TrackRef> & input);
    bool select(reco::TrackCollection & nonTauTracks, vVVTrackRef & threeProngCombis);
	bool removeDuplicateTriplets(const std::vector<reco::TrackRef> & duplicateTracks, vVVTrackRef & threeProngCombis, vVVTrackRef::iterator & candidates, vVTrackRef::iterator & triplets);		
	
	edm::Event * iEvent_;
    edm::ParameterSet iConfig_;
	edm::InputTag inputCollectionTag_, trackCollectionTag_;
	unsigned int cnt_, cntFound_, minTau_;
    
    template <typename T> double getInvariantMass(const T & tracks, const double mass = 0.140){//if second argument empty default pion is supposed
		double SumPx = 0;
		double SumPy = 0;
		double SumPz = 0;
		double SumE = 0;
		
		for(unsigned int i=0; i<tracks.size(); i++){
			SumPx += tracks[i]->px();
			SumPy += tracks[i]->py();
			SumPz += tracks[i]->pz();
			SumE += sqrt(pow(tracks[i]->p(),2)+pow(mass,2));
		}
		double invmass = sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2));
		return invmass;
	}
	template <typename T> static bool cmpPt(const T & a, const T & b){
		return a->pt() > b->pt();
	}
};
