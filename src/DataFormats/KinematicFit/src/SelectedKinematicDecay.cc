#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"



SelectedKinematicDecay::SelectedKinematicDecay() {

}
SelectedKinematicDecay::SelectedKinematicDecay(SelectedKinematicParticleCollection particles)
{
    particles_ = particles;
	signalPFChargedHadrCands_ = -1;
	signalPFNeutrHadrCands_ = -1;
	discriminators_.clear();
}
SelectedKinematicDecay::SelectedKinematicDecay(SelectedKinematicParticleCollection particles, const int & signalPFChargedHadrCands, const int & signalPFNeutrHadrCands, const std::map<std::string, bool> & discriminators)
{
    particles_ = particles;
	signalPFChargedHadrCands_ = signalPFChargedHadrCands;
	signalPFNeutrHadrCands_   = signalPFNeutrHadrCands;
	discriminators_ = discriminators;
}

const SelectedKinematicParticle* SelectedKinematicDecay::topParticle() const
{
    return &(particles_.front());
}
void SelectedKinematicDecay::particles(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		par.push_back(&(*iter));
    }
}
void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
}
void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 0 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
int SelectedKinematicDecay::signalPFChargedHadrCands(){
	return signalPFChargedHadrCands_;
}
int SelectedKinematicDecay::signalPFNeutrHadrCands(){
	return signalPFNeutrHadrCands_;
}
std::map<std::string, bool> SelectedKinematicDecay::discriminators() const
{
	return discriminators_;
}
void SelectedKinematicDecay::modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par)
{
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
