#ifndef ACPFTauConverter_h
#define ACPFTauConverter_h
/**
 PFTau helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACDataFormats/ACGeneral/interface/ACPFTau.h"
#include "ACFrameworkModules/Common/interface/ConversionLog.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

class ACPFTauConverter : public ACPFTau {
public:
    virtual std::string classname() const { return "ACPFTauConverter"; }
    /// constructor from PFTau to initiate a ACPFTau
    ACPFTauConverter(const edm::Event& evt, const reco::PFTauRef& tau, const std::vector<edm::InputTag> & discriminators, const PFJetConversionLog * conversionLogPFJet);
    virtual ~ACPFTauConverter() {};
};

#endif
