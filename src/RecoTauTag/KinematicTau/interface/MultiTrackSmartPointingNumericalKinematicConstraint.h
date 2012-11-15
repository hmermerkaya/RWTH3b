#ifndef MultiTrackSmartPointingNumericalKinematicConstraint_H
#define MultiTrackSmartPointingNumericalKinematicConstraint_H

#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"

class MultiTrackSmartPointingNumericalKinematicConstraint : public MultiTrackNumericalKinematicConstraint {
 public:
  MultiTrackSmartPointingNumericalKinematicConstraint(GlobalPoint& ref, double weight=1.0):MultiTrackNumericalKinematicConstraint(weight),refPoint(ref){};

  virtual int numberOfEquations() const {return 2;}

  virtual MultiTrackSmartPointingNumericalKinematicConstraint *clone() const
  {return new MultiTrackSmartPointingNumericalKinematicConstraint(*this);}

    
 protected:
  virtual AlgebraicVector  Compute(const std::vector<KinematicState> states,const GlobalPoint& point, const std::vector<double> &epsilonPar, const std::vector<double> &epsilonPos) const;
  
 private:
  GlobalPoint refPoint;

};
#endif
