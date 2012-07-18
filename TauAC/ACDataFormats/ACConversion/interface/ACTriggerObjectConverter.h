#ifndef ACTriggerObjectConverter_h
#define ACTriggerObjectConverter_h
/**
 Helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "../../ACGeneral/interface/ACTriggerObject.h"

class ACTriggerObjectConverter : public ACTriggerObject {
public:
    virtual std::string classname() const { return "ACTriggerObjectConverter"; }
    /// constructor from trigger::TriggerObject
    ACTriggerObjectConverter(const trigger::TriggerObject & triggerObject);
    virtual ~ACTriggerObjectConverter() {};//do we need a destructor here?
};

#endif
