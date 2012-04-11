#ifndef ACFittedThreeProngDecayConverter_h
#define ACFittedThreeProngDecayConverter_h
/**
 Decay helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "ACFrameworkModules/Common/interface/MCTruthMatching.h"
#include "ACFrameworkModules/Common/interface/ConversionLog.h"
#include "ACDataFormats/ACGeneral/interface/ACFitParticle.h"
#include "ACDataFormats/ACGeneral/interface/ACFittedThreeProngDecay.h"
#include "ACDataFormats/ACGeneral/interface/ACGenDecay.h"
#include "ACFitParticleConverter.h"
#include "DataFormats/Math/interface/deltaR.h"

class ACFittedThreeProngDecayConverter : public ACFittedThreeProngDecay {
public:
    virtual std::string classname() const { return "ACFittedThreeProngDecayConverter"; }
    /// constructor from SelectedKinematicDecay to initiate a ACFittedThreeProngDecay
    ACFittedThreeProngDecayConverter(const edm::Event& evt, const SelectedKinematicDecay & decay, MCTruthMatching * kinematicParticleMatching, std::vector<ACGenDecay *> * genTauDecays, std::vector<ACFitParticle *> * kinematicParticles,  PFTauConversionLog * conversionLogPFTau);
    virtual ~ACFittedThreeProngDecayConverter() {};

};

#endif
