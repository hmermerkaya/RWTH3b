#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"



SelectedKinematicDecay::SelectedKinematicDecay() {

}
SelectedKinematicDecay::SelectedKinematicDecay(SelectedKinematicParticleCollection particles)
{
    particles_ = particles;
}

SelectedKinematicParticle* SelectedKinematicDecay::topParticle()
{
    return &(particles_.front());
}
void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle* > & par)
{
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
}
void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle* > & par)
{
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
//std::vector< SelectedKinematicParticle* > SelectedKinematicDecay::chargedDaughters()
//{
//    std::vector< SelectedKinematicParticle* > tmpVec;
//    for ( unsigned int i=0; i!=particles_.size(); i++ ) {
//        if ( std::abs(particles_[i].charge()) == 1 ) {
//			if(i != 0) tmpVec.push_back(&(particles_[i]));//skip mother
//        }
//    }
//    return tmpVec;
//}


void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle* > & par)
{
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 0 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
