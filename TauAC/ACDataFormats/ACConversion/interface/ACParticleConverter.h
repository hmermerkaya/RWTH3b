#ifndef ACParticleConverter_h
#define ACParticleConverter_h
/**
 Particle helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class ACParticleConverter : public ACParticle {
public:
    virtual std::string classname() const { return "ACParticleConverter"; }
    /// constructor from SelectedKinematicParticle to initiate a ACParticle
    ACParticleConverter(const reco::Muon& particle, ACGenParticleRef genRef = ACGenParticleRef());
    ACParticleConverter(const reco::GsfElectron& particle, ACGenParticleRef genRef = ACGenParticleRef());
    virtual ~ACParticleConverter() {};
};

#endif
