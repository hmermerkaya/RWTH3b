#ifndef ACGenDecay_h
#define ACGenDecay_h
/**
 Class to store decay information. The template argument expects either T=ACParticle or T=ACFitParticle

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACDecay.h"
#include "ACGenParticle.h"

class ACGenDecay : public ACDecayBase<ACGenParticle> {
public:
    virtual std::string classname() const { return "ACGenDecay"; }
    ACGenDecay();
    ACGenDecay(std::vector<ACGenParticleRef> particles, std::string name, ACGenDecayRef genDecayRef = ACGenDecayRef());
    virtual ~ACGenDecay() {};

    ClassDef(ACGenDecay,1);

protected:
    /// no new member for now, but need a class to be able to store to file
};

typedef ACGenDecay* pACGenDecay;

#endif
