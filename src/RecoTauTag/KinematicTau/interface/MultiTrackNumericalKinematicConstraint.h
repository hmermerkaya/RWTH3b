#ifndef MultiTrackNumericalKinematicConstraint_H
#define MultiTrackNumericalKinematicConstraint_H

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"

class MultiTrackNumericalKinematicConstraint {
 public:
  enum Position{pos_x=0,pos_y,pos_z,nposdim};
  enum Parameters{par_vx=0,par_vy,par_vz,par_px,par_py,par_pz,par_m,npardim};

  MultiTrackNumericalKinematicConstraint(bool debugflag,double weight=1.0);
  virtual ~MultiTrackNumericalKinematicConstraint(){};

  virtual MultiTrackNumericalKinematicConstraint* clone() const=0;

  virtual double GetWeight()const{return weight_;}

  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint)=0;
  bool  ApplyLagrangianConstraints(double &chi_2, double &delta);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint)=0;
  virtual int numberOfEquations()=0;
  virtual bool isConverged(double chi_2, double delta);

 protected:
  virtual TVectorD Value(TVectorD &v)=0;
  
  TVectorD par_first;
  TMatrixTSym<double> cov_first;
  TVectorD par;
  TMatrixTSym<double> cov;
  double chi2;
  TVectorD maxStep;
  TVectorD w; // weight for constraints
  const MagneticField* field;
  bool debug;

 private:
  TMatrixT<double> Derivative();
  TVectorT<double> convertToVector(TMatrixT<double> M);
  TMatrixT<double> convertToMatrix(TVectorT<double> V);

  double epsilon_, weight_;
  
};
#endif
