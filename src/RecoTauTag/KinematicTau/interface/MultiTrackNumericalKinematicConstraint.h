#ifndef MultiTrackNumericalKinematicConstraint_H
#define MultiTrackNumericalKinematicConstraint_H

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

class MultiTrackNumericalKinematicConstraint {
 public:
  enum ParameterSize{npardim=7,nposdim=3};

  MultiTrackNumericalKinematicConstraint(double weight=1.0);
  virtual ~MultiTrackNumericalKinematicConstraint(){};

  virtual MultiTrackNumericalKinematicConstraint* clone() const=0;

  virtual double GetWeight()const{return weight_;}

  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint)=0;
  bool  ApplyLagrangianConstraints(double &chi_2, double &delta);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint)=0;
  virtual int numberOfEquations()=0;

 protected:
  virtual AlgebraicVector Value(AlgebraicVector &v)=0;
  
  AlgebraicVector par_first;
  AlgebraicMatrix cov_first;
  AlgebraicVector par;
  AlgebraicMatrix cov;
  AlgebraicVector chi2;
  const MagneticField* field;

 private:
  AlgebraicMatrix Derivative();

  double epsilon_, weight_;
  
};
#endif
