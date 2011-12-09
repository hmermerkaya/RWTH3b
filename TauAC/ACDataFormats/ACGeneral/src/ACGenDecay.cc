#include "../interface/ACGenDecay.h"

ClassImp(ACGenDecay)

ACGenDecay::ACGenDecay():
ACDecayBase<ACGenParticle>() {
}
ACGenDecay::ACGenDecay(std::vector<ACGenParticleRef> particles, std::string name, ACGenDecayRef genDecayRef):
ACDecayBase<ACGenParticle>(particles, name, genDecayRef) {
}
