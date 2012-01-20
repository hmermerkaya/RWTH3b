#ifndef ACGenParticleConverter_h
#define ACGenParticleConverter_h
/**
 Particle helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACGenParticle.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/Candidate/interface/Candidate.h>

class ACGenParticleConverter : public ACGenParticle {
public:
    virtual std::string classname() const { return "ACGenParticleConverter"; }
    /// constructor from SelectedKinematicParticle to initiate a ACGenParticle
    ACGenParticleConverter(const reco::GenParticle& particle);
    virtual ~ACGenParticleConverter() {};

    /// store the reference to the generator tau decay that contains the current particle
    void setGenDecayRef(const ACGenDecayRef & genDecayRef);
    /// store the reference to the generated mother particle (if stored)
    void setMother(const ACGenParticleRef & motherRef);
    /// store the reference to a generated daughter particle (if stored)
    void addDaughter(const ACGenParticleRef & daughterRef);
};

#endif
