#ifndef MultiTrackNumericalKinematicConstraint_H
#define MultiTrackNumericalKinematicConstraint_H

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"

class MultiTrackNumericalKinematicConstraint : public MultiTrackKinematicConstraint{

 public:
  MultiTrackNumericalKinematicConstraint(double weight=1.0);

  virtual AlgebraicVector  value(const std::vector<KinematicState> states,const GlobalPoint& point) const;

  virtual AlgebraicMatrix parametersDerivative(const std::vector<KinematicState> states,const GlobalPoint& point) const;

  virtual AlgebraicMatrix positionDerivative(const std::vector<KinematicState> states,const GlobalPoint& point) const;
  
  virtual double GetWeight()const{return weight_;}

 protected:
  virtual AlgebraicVector Compute(const std::vector<KinematicState> states,const GlobalPoint& point, const std::vector<double> &epsilonPar, const std::vector<double> &epsilonPos) const=0;
  unsigned int npardim,nposdim;

 private:
  virtual AlgebraicVector  ComputeDefault(const std::vector<KinematicState> states,const GlobalPoint& point) const;
  double distance_epsilon,momentum_epsilon;
  double weight_;
};
#endif
