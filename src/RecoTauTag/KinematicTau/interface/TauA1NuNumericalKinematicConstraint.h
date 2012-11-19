#ifndef TauA1NuNumericalKinematicConstraint_H
#define TauA1NuNumericalKinematicConstraint_H

#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "TVector3.h"

class TauA1NuNumericalKinematicConstraint : public MultiTrackNumericalKinematicConstraint{
 public:
  TauA1NuNumericalKinematicConstraint(TVector3 pv,double weight=1.0):
    MultiTrackNumericalKinematicConstraint(weight)
    {
      pv_inital=pv;
    }
  virtual ~TauA1NuNumericalKinematicConstraint(){}

  virtual TauA1NuNumericalKinematicConstraint * clone() const {return new TauA1NuNumericalKinematicConstraint(*this);}
  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint);
  virtual int numberOfEquations(){return 3;}

 protected:
  virtual AlgebraicVector Value(AlgebraicVector &v);
    
 private:
  TVector3 pv_inital;
};
#endif
