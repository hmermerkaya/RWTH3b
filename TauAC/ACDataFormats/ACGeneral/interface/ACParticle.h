#ifndef ACParticle_h
#define ACParticle_h
/**
 Class to store particle information. Adds additional members (like covariance matrices) to ACGenParticles.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACGenParticle.h"
#include "ACVertex.h"
#include "ACRef.h"

class ACParticle : public ACGenParticle {
public:
    virtual std::string classname() const { return "ACParticle"; }
    ACParticle();
    ACParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex = TVector3(-1000,-1000,-1000), const TMatrixDSym & matrix = TMatrixDSym(7), ACGenParticleRef genRef = ACGenParticleRef(), ACGenDecayRef genDecayRef = ACGenDecayRef());
    virtual ~ACParticle() {};

    const TMatrixDSym & matrix() const;
    const ACGenParticleRef & genRef() const;

    /// vertex including the corresponding covariance
    const ACVertex vtx() const;//fixme: find better name (conflict with 3d point in ACCandidate)
    /// vertex part of the covariance matrix
    const TMatrixDSym vertexError() const;
    /// corresponding diagonal element of the covariance matrix
    double svx() const;
    /// corresponding diagonal element of the covariance matrix
    double svy() const;
    /// corresponding diagonal element of the covariance matrix
    double svz() const;
    /// corresponding diagonal element of the covariance matrix
    double spx() const;
    /// corresponding diagonal element of the covariance matrix
    double spy() const;
    /// corresponding diagonal element of the covariance matrix
    double spz() const;
    /// corresponding diagonal element of the covariance matrix
    double sm() const;
    /// error on the transversal momentum
    double spt() const;
    /// error on the momentum magnitude
    double sp() const;

    ClassDef(ACParticle,1);

protected:
    /// covariance matrix 7x7, according to the KinematicFit package (x, y, z, px, py, pz, m)
    TMatrixDSym matrix_;
    /// reference to the matched generated particle
    ACGenParticleRef genRef_;
};

typedef ACParticle* pACParticle;

#endif
