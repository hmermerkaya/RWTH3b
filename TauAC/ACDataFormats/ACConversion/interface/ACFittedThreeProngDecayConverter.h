#ifndef ACFittedThreeProngDecayConverter_h
#define ACFittedThreeProngDecayConverter_h
/**
 Decay helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "ACFrameworkModules/Common/interface/MCTruthMatching.h"
#include "ACFrameworkModules/Common/interface/PFTauMatching.h"
#include "../../ACGeneral/interface/ACFitParticle.h"
#include "../../ACGeneral/interface/ACFittedThreeProngDecay.h"
#include "ACFitParticleConverter.h"

class ACFittedThreeProngDecayConverter : public ACFittedThreeProngDecay {
public:
    virtual std::string classname() const { return "ACFittedThreeProngDecayConverter"; }
    /// constructor from SelectedKinematicDecay to initiate a ACFittedThreeProngDecay
    ACFittedThreeProngDecayConverter(const edm::Event& evt, const SelectedKinematicDecay & decay, MCTruthMatching * kinematicParticleMatching, std::vector<ACFitParticle *> * kinematicParticles,  PFTauMatching * pfTauMatching);
    virtual ~ACFittedThreeProngDecayConverter() {};

};

#endif
