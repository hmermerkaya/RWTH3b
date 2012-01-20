#ifndef ACGenParticle_h
#define ACGenParticle_h
/**
 Base class to store generator-particle information. Mother/daughter relations are maintained.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TMath.h>
#include <TVector3.h>
#include <TMatrixDSym.h>

#include "ACCandidate.h"
#include "ACRef.h"

class ACGenParticle : public ACCandidate {
public:
    virtual std::string classname() const { return "ACGenParticle"; }
    ACGenParticle();
    ACGenParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex = TVector3(-1000,-1000,-1000), const ACGenDecayRef & genDecayRef = ACGenDecayRef(), const ACGenParticleRef & mother = ACGenParticleRef(), const std::vector<ACGenParticleRef> & daughters = std::vector<ACGenParticleRef>());
    virtual ~ACGenParticle() {};

    int status() const;
    int pdgId() const;
    const TVector3 & vertex() const;
    const ACGenDecayRef & genDecayRef() const;
    const ACGenParticleRef & mother() const;
    const std::vector<ACGenParticleRef> & daughters() const;

    /// x coordinate of vertex position
    double vx() const;
    /// y coordinate of vertex position
    double vy() const;
    /// z coordinate of vertex position
    double vz() const;

    /// q/p of the Helix parametrization
    double qoverp() const;
    /// lambda of the Helix parametrization
    double lambda() const;
    /// dxy of the Helix parametrization
    double dxy() const;
    /// dsz of the Helix parametrization
    double dsz() const;

    ClassDef(ACGenParticle,1);

protected:
    /// status code
    int status_;
    /// PDG identifier (assumption in case of reco)
    int pdgId_;
    /// vertex position
    TVector3 vertex_;
    /// reference to the generator tau decay that contains this particle (invalid ref if not contained in any)
    ACGenDecayRef genDecayRef_;
    /// reference to the mother particle (if any)
    ACGenParticleRef mother_;
    /// reference to daughter particles (if any)
    std::vector<ACGenParticleRef> daughters_;
};

typedef ACGenParticle* pACGenParticle;

#endif
