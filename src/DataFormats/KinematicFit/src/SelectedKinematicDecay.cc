#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"



SelectedKinematicDecay::SelectedKinematicDecay() {
    particles_.clear();
    iterations_ = -1;
    maxiterations_ = -1;
    csum_ = -1.;
    mincsum_ = -1.;
    chi2_ = -1.; 
	constraints_ = 0.;
	ndf_ = 0.;
    PFTauRef_ = reco::PFTauRef();
    discriminators_.clear();
    //    primaryVertexRef_ = reco::VertexRef();
    
    vtxSignPVRotSV_ = -1.;
    vtxSignPVRotPVRed_ = -1.;
    a1Mass_ = -1.;
    energyTFraction_ = -1.;
}
SelectedKinematicDecay::SelectedKinematicDecay(const SelectedKinematicParticleCollection & particles, const int iterations, const int maxiterations, const float csum, const float mincsum, const int constraints, const int ndf, const float chi2, const reco::PFTauRef & tauRef, const std::map<std::string, bool> & discriminators/*, const reco::VertexRef & primaryVertexRef*/) {
    particles_ = particles;
    iterations_ = iterations;
    maxiterations_ = maxiterations;
    csum_ = csum;
    mincsum_ = mincsum;
    constraints_ = constraints;
	ndf_ = ndf;
    chi2_ = chi2;
    
    PFTauRef_ = tauRef;
	discriminators_ = discriminators;
    //    primaryVertexRef_ = primaryVertexRef;
}

const SelectedKinematicParticle * SelectedKinematicDecay::topParticle() const {
    return &(particles_.front());
}
void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle const * > & par) const {
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
}
void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const {
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const {
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 0 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}

void SelectedKinematicDecay::modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par) {
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}

void SelectedKinematicDecay::setMissingQualityCriteria(const double & vtxSignPVRotSV, const double & vtxSignPVRotPVRed, const double & a1Mass, const double & energyTFraction) {
    vtxSignPVRotSV_ = vtxSignPVRotSV;
    vtxSignPVRotPVRed_ = vtxSignPVRotPVRed;
    
    //WARNING!!!
    //from now one we assume a tau decay into three pions and neutrino
    //other channels need their own discriminators
    //!!!
    a1Mass_ = a1Mass;
    energyTFraction_ = energyTFraction;
}
const int SelectedKinematicDecay::iterations() const {
    return iterations_;
}
const int SelectedKinematicDecay::maxiterations() const {
    return maxiterations_;
}
const float SelectedKinematicDecay::chi2() const {
    return chi2_;
}
const float SelectedKinematicDecay::constraints() const {
    return constraints_;
}
const float SelectedKinematicDecay::ndf() const {
    return ndf_;
}
const float SelectedKinematicDecay::csum() const {
    return csum_;
}
const float SelectedKinematicDecay::mincsum() const {
    return mincsum_;
}

const double SelectedKinematicDecay::chi2prob() const {
    ChiSquared chiSquared(chi2(), ndf());
    return chiSquared.probability();
}
const int SelectedKinematicDecay::sgnlConeTrkSize() const {
    if (PFTauRef().isAvailable()) {
        return PFTauRef()->signalPFChargedHadrCands().size();
    }
    printf("SelectedKinematicDecay::sgnlConeTrkSize:ERROR! reference to PFTau invalid. Return dummy value.\n");
    return 0;
}
