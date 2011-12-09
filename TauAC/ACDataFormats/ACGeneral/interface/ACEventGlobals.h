#ifndef ACEventGlobals_h
#define ACEventGlobals_h
/**
 Class containing global event data.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TVector3.h>

#include "ACCommon.h"

class ACEventGlobals : public ACCommon {
public:
    virtual std::string classname() const { return "ACEventGlobals"; }
    ACEventGlobals(const TVector3& pfMET, const TVector3& tcMET, double pfSumET, double tcSumET);
    ACEventGlobals();
    virtual ~ACEventGlobals() {};

    const TVector3& pfMET() const { return pfMET_; }
    const TVector3& tcMET() const { return tcMET_; }
    double pfSumET() const { return pfSumET_; }
    double tcSumET() const { return tcSumET_; }
    double pfMETSignificance() const;
    double tcMETSignificance() const;
    double pfMETTopology() const;
    double tcMETTopology() const;

    ClassDef(ACEventGlobals,1);

protected:
    /// missing transverse energy
    TVector3 pfMET_, tcMET_;
    /// total transverse energy
    double pfSumET_, tcSumET_;
};

typedef ACEventGlobals* pACEventGlobals;

#endif
