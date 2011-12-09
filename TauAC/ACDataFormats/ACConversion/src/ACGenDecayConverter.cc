#include "../interface/ACGenDecayConverter.h"

ACGenDecayConverter::ACGenDecayConverter(const ACGenParticleRef & topParticle) {
    /// genDecayRef_ won't be set
    particles_.clear();
    particles_.push_back(topParticle);
    if (topParticle->charge() == -1) {
        name_ = std::string("Tau-");
    } else {
        name_ = std::string("Tau+");
    }
}

void ACGenDecayConverter::addDaughter(const ACGenParticleRef & daughter) {
    particles_.push_back(daughter);
}
