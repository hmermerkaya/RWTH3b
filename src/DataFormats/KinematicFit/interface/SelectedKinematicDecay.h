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
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles, const int & signalPFChargedHadrCands, const int & signalPFNeutrHadrCands);

    const SelectedKinematicParticle* topParticle() const;
    void particles(std::vector< SelectedKinematicParticle const * > & par) const;
    void daughters(std::vector< SelectedKinematicParticle const * > & par) const;
    void chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
    void neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
	int signalPFChargedHadrCands();
	int signalPFNeutrHadrCands();
	/**
	 DO NOT USE after reading from event stream!
	 */
	void modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par);
 	
private:
    SelectedKinematicParticleCollection particles_;
	/**
	 size of tracks in signal cone of PFTau candidate
	 */
	int signalPFChargedHadrCands_;
	int signalPFNeutrHadrCands_;
};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
