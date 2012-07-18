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
    ACEventGlobals(const TVector3& pfMET, const TVector3& pfType1CorrectedMET, const TVector3& tcMET, double pfSumET, double pfType1CorrectedSumET, double tcSumET);
    ACEventGlobals();
    virtual ~ACEventGlobals() {};

    const TVector3& pfMET() const { return pfMET_; }
    const TVector3& pfType1CorrectedMET() const { return pfType1CorrectedMET_; }
    const TVector3& tcMET() const { return tcMET_; }
    double pfSumET() const { return pfSumET_; }
    double pfType1CorrectedSumET() const { return pfType1CorrectedSumET_; }
    double tcSumET() const { return tcSumET_; }
    double pfMETSignificance() const;
    double tcMETSignificance() const;
    double pfMETTopology() const;
    double tcMETTopology() const;

    ClassDef(ACEventGlobals,1);

protected:
    /// missing transverse energy from PFJets
    TVector3 pfMET_;
    /// missing transverse energy from PFJets type-1 corrected
    TVector3 pfType1CorrectedMET_;
    /// missing transverse energy from CaloJets
    TVector3 tcMET_;
    /// total transverse energy from PFJets
    double pfSumET_;
    /// total transverse energy from PFJets type-1 corrected
    double pfType1CorrectedSumET_;
    /// total transverse energy from CaloJets
    double tcSumET_;
};

typedef ACEventGlobals* pACEventGlobals;

#endif
