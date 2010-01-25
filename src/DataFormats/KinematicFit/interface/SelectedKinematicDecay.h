#ifndef DataFormats_KinematicFit_SelectedKinematicDecay_h
#define DataFormats_KinematicFit_SelectedKinematicDecay_h
// -*- C++ -*-
//
// Package:    SelectedKinematicDecay
// Class:      SelectedKinematicDecay
// 
/**
 
 Description: saves SelectedKinematicParticle of one tau decay
 
 Implementation:
 <Notes on implementation>
 */
//
//

#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"



class SelectedKinematicDecay {
public:
	SelectedKinematicDecay();
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles);

    SelectedKinematicParticle* topParticle();
    void daughters(std::vector< SelectedKinematicParticle* > & par);
    void chargedDaughters(std::vector< SelectedKinematicParticle* > & par);
    void neutralDaughters(std::vector< SelectedKinematicParticle* > & par);
 	
private:
    SelectedKinematicParticleCollection particles_;
};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
