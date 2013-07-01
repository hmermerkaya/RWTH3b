#ifndef MultiTrackMassNumericalKinematicConstraint_H
#define MultiTrackMassNumericalKinematicConstraint_H

#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"

/**
 * Constraint to force some of the particles in the fit to have a certain invariant mass.
 */

class MultiTrackMassNumericalKinematicConstraint : public MultiTrackNumericalKinematicConstraint {
 public:
  MultiTrackMassNumericalKinematicConstraint(const ParticleMass& theMass, const unsigned int nbrParticles, double weight=1.0):MultiTrackNumericalKinematicConstraint(weight), mass(theMass), nPart(nbrParticles){};

 virtual int numberOfEquations() const {return 1;}

 virtual MultiTrackMassNumericalKinematicConstraint *clone() const
 {return new MultiTrackMassNumericalKinematicConstraint(*this);}

    
 protected:
  virtual AlgebraicVector  Compute(const std::vector<KinematicState> states,const GlobalPoint& point, const std::vector<double> &epsilonPar, const std::vector<double> &epsilonPos) const;
  
 private:
  const ParticleMass mass;
  const unsigned int nPart;

};
#endif
