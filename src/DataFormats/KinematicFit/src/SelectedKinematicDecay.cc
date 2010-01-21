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
std::vector< SelectedKinematicParticle* > SelectedKinematicDecay::daughters()
{
    std::vector< SelectedKinematicParticle* > tmpVec;
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        tmpVec.push_back(&(*iter));
    }
    return tmpVec; 
}
std::vector< SelectedKinematicParticle* > SelectedKinematicDecay::chargedDaughters()
{
    std::vector< SelectedKinematicParticle* > tmpVec;
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
            tmpVec.push_back(&(*iter));
        }
    }
    return tmpVec;
}
std::vector< SelectedKinematicParticle* > SelectedKinematicDecay::neutralDaughters()
{
    std::vector< SelectedKinematicParticle* > tmpVec;
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 0 ) {
            tmpVec.push_back(&(*iter));
        }
    }
    return tmpVec;
}
