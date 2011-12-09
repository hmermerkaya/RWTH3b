#ifndef ACCandidate_h
#define ACCandidate_h
/**
 Class to store generic candidate objects.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TLorentzVector.h>

#include "ACCommon.h"

class ACCandidate : public ACCommon {
public:
    virtual std::string classname() const { return "ACCandidate"; }
    ACCandidate();
    ACCandidate(const TLorentzVector& p4, float charge_ = 0.0);
    virtual ~ACCandidate() {};

    float charge() const;
    const TLorentzVector & p4() const;

    double p() const;
    /// energy
    double energy() const;
    /// transverse energy
    double et() const;
    /// mass
    double mass() const;
    /// transverse mass
    double mt() const;
    /// x coordinate of momentum vector
    double px() const;
    /// y coordinate of momentum vector
    double py() const;
    /// z coordinate of momentum vector
    double pz() const;
    /// transverse momentum
    double pt() const;
    /// momentum azimuthal angle
    double phi() const;
    /// momentum polar angle
    double theta() const;
    /// momentum pseudorapidity
    double eta() const;

    ClassDef(ACCandidate,1);

protected:
    /// four-momentum Lorentz vector
    TLorentzVector p4_;
    /// charge
    float charge_;
};

typedef ACCandidate* pACCandidate;

#endif
