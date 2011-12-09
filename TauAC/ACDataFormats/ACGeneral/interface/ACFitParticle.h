#ifndef ACFitParticle_h
#define ACFitParticle_h
/**
 Derived class to store additional particle properties occuring in the KinematicFit (like an initial state).

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACParticle.h"
#include "ACRef.h"

class ACFitParticle : public ACParticle {
public:
    const std::string classname() { return "ACFitParticle"; }
    ACFitParticle();
    ACFitParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex, const TMatrixDSym & matrix, const TLorentzVector & initialP4, const TVector3 & initialVertex, const TMatrixDSym & initialMatrix, ACGenParticleRef genRef = ACGenParticleRef(), ACGenDecayRef genDecayRef = ACGenDecayRef());
    virtual ~ACFitParticle() {};

    const TLorentzVector & initialP4() const;
    const TVector3 & initialVertex() const;
    const TMatrixDSym & initialMatrix() const;

    /// vertex part of the initial covariance matrix
    const TMatrixDSym vertexError() const;
    /// corresponding diagonal element of the initial covariance matrix
    double svxInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double svyInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double svzInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double spxInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double spyInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double spzInitial() const;
    /// corresponding diagonal element of the initial covariance matrix
    double smInitial() const;

    /// q/p of the Helix parametrization
    double qoverpInitial() const;
    /// lambda of the Helix parametrization
    double lambdaInitial() const;
    /// dxy of the Helix parametrization
    double dxyInitial() const;
    /// dsz of the Helix parametrization
    double dszInitial() const;

    ClassDef(ACFitParticle,1);

protected:
    /// initial four-momentum Lorentz vector before fit
    TLorentzVector initialP4_;
    /// initial vertex position before fit
    TVector3 initialVertex_;
    /// covariance matrix 7x7 of the initial state, according to the KinematicFit package (x, y, z, px, py, pz, m)
    TMatrixDSym initialMatrix_;
};

typedef ACFitParticle* pACFitParticle;

#endif
