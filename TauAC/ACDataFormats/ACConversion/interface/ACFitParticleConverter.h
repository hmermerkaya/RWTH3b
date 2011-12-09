#ifndef ACFitParticleConverter_h
#define ACFitParticleConverter_h
/**
 Particle helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"
#include "../../ACGeneral/interface/ACFitParticle.h"

class ACFitParticleConverter : public ACFitParticle {
public:
    virtual std::string classname() const { return "ACFitParticleConverter"; }
    /// constructor from SelectedKinematicParticle to initiate a ACFittedThreeProngDecay
    ACFitParticleConverter(const SelectedKinematicParticle & selPart, int decayCharge, const ACGenParticleRef * genRef);
    virtual ~ACFitParticleConverter() {};//do we need a destructor here?

private:
    int getPdgId(const std::string name, const int charge) const;
};

#endif
