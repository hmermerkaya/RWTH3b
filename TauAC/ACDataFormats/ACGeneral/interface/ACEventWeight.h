#ifndef ACEventWeight_h
#define ACEventWeight_h
/**
 Class containing the event weight.

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

    double weight() const { return weight_; }
    
    ClassDef(ACEventWeight,1);

protected:
    /// event weight
    double weight_;
};

typedef ACEventWeight* pACEventWeight;

#endif
