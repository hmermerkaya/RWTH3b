#ifndef ACJetConverter_h
#define ACJetConverter_h
/**
 Jet helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACJet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

class ACJetConverter : public ACJet {
public:
    virtual std::string classname() const { return "ACJetConverter"; }
    /// constructor from PFJet to initiate a ACJet
    ACJetConverter(const reco::Jet& jet);
    ACJetConverter(const reco::PFJet& jet);
    virtual ~ACJetConverter() {};
};

#endif
