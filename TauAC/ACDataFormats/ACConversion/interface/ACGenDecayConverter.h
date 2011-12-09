#ifndef ACGenDecayConverter_h
#define ACGenDecayConverter_h
/**
 Decay helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACGenDecay.h"
#include "ACFitParticleConverter.h"

class ACGenDecayConverter : public ACGenDecay {
public:
    virtual std::string classname() const { return "ACGenDecayConverter"; }
    /// constructor from the top particle of the decay
    ACGenDecayConverter(const ACGenParticleRef & topParticle);
    virtual ~ACGenDecayConverter() {}; //do we need a destructor here?

    /// add daughters to the decay
    void addDaughter(const ACGenParticleRef & daughter);

};

#endif
