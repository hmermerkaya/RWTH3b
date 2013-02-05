#ifndef ACEventGlobals_h
#define ACEventGlobals_h
/**
 Class containing global event data.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TVector3.h>

#include "ACCommon.h"
#include "ACMET.h"

class ACEventGlobals : public ACCommon {
public:
    virtual std::string classname() const { return "ACEventGlobals"; }
    ACEventGlobals(){}
    ACEventGlobals(const ACMET & pfMET, const ACMET & pfType1CorrectedMET, const ACMET & tcMET);
    virtual ~ACEventGlobals() {};

    const ACMET & pfMET() const { return pfMET_; }
    const ACMET & pfType1CorrectedMET() const { return pfType1CorrectedMET_; }
    const ACMET & tcMET() const { return tcMET_; }

    ClassDef(ACEventGlobals,1);

protected:
    /// missing transverse energy from PFJets
    ACMET pfMET_;
    /// missing transverse energy from PFJets type-1 corrected
    ACMET pfType1CorrectedMET_;
    /// missing transverse energy from CaloJets
    ACMET tcMET_;
};

typedef ACEventGlobals* pACEventGlobals;

#endif
