#ifndef ACTrack_h
#define ACTrack_h
/**
 Base class to store (lightweight) track information (similar to CMSSW's reco::Track).
 Move to a real helix parametrization?

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TVector3.h>

#include "ACCandidate.h"

class ACTrack : public ACCandidate {
public:
    virtual std::string classname() const { return "ACTrack"; }
    ACTrack();
    ACTrack(const TLorentzVector & p4, float charge, const std::string & name, float chi2, float ndof, const TVector3 & vertex = TVector3(-1000,-1000,-1000));
    virtual ~ACTrack() {};

    const TVector3 & vertex() const;
    float chi2() const { return chi2_; }
    float ndof() const { return ndof_; }
    /// x coordinate of vertex position
    double vx() const;
    /// y coordinate of vertex position
    double vy() const;
    /// z coordinate of vertex position
    double vz() const;

    ClassDef(ACTrack,1);

protected:
    /// chi2 of the track fit
    float chi2_;
    /// ndof of the track fit
    float ndof_;
    /// vertex position
    TVector3 vertex_;
};

typedef ACTrack* pACTrack;

#endif
