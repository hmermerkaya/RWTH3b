#ifndef ACEventWeight_h
#define ACEventWeight_h
/**
 Class containing event weights.

 @author Lars Perchalla & Philip Sauerland
 @date 2012
 */

#include "ACCommon.h"

class ACEventWeight : public ACCommon {
public:
    virtual std::string classname() const { return "ACEventWeight"; }
    ACEventWeight(double weight);
    ACEventWeight();
    virtual ~ACEventWeight() {};

    const double & weight() const { return weight_; }
    
    ClassDef(ACEventWeight,1);

protected:
    /// event weight derived from pileup reweighting
    double weight_;
    //fixme: rename member once additional weights are needed
};

typedef ACEventWeight* pACEventWeight;

#endif
