#ifndef DataFormats_KinematicFit_SelectedKinematicDecay_h
#define DataFormats_KinematicFit_SelectedKinematicDecay_h
// -*- C++ -*-
//
// Package:    SelectedKinematicDecay
// Class:      SelectedKinematicDecay
// 
/**
 
 Description: stores SelectedKinematicParticle of one tau decay
 
 Implementation:
 <Notes on implementation>
 */
//
//
// Original Author:  Lars Perchalla, Philip Sauerland


#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"



class SelectedKinematicDecay {
public:
	SelectedKinematicDecay();
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles);
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles, const int & signalPFChargedHadrCands);

    SelectedKinematicParticle* topParticle();
    void particles(std::vector< SelectedKinematicParticle* > & par);
    void daughters(std::vector< SelectedKinematicParticle* > & par);
    void chargedDaughters(std::vector< SelectedKinematicParticle* > & par);
    void neutralDaughters(std::vector< SelectedKinematicParticle* > & par);
	int signalPFChargedHadrCands();
 	
private:
    SelectedKinematicParticleCollection particles_;
	/**
	 size of tracks in signal cone of PFTau candidate
	 */
	int signalPFChargedHadrCands_;
};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
