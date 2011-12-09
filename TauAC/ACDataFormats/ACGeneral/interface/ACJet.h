#ifndef ACJet_h
#define ACJet_h
/**
 Class to store jets.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCandidate.h"

class ACJet : public ACCandidate {
public:
    virtual std::string classname() const { return "ACJet"; }
    ACJet();
    ACJet(const TLorentzVector& p4, int multiplicity, int charge);
    virtual ~ACJet() {};

    int multiplicity() const { return multiplicity_; }

    ClassDef(ACJet,1);

protected:
    /// jet multiplicity
    int multiplicity_;
};

typedef ACJet* pACJet;

#endif
