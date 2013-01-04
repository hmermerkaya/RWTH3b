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
  enum ConvergeProc{ConstraintMin=0,Chi2Min,Chi2AndConstaintMin};

  MultiTrackNumericalKinematicConstraint(bool debugflag,double weight);
  virtual ~MultiTrackNumericalKinematicConstraint(){};

  virtual MultiTrackNumericalKinematicConstraint* clone() const=0;

  virtual void   SetWeight(double weight){weight_=weight;}
  virtual double GetWeight()const{return weight_;}
  virtual void SetMaxDelta(double MaxDelta){MaxDelta_=MaxDelta;}

  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint)=0;
  bool  ApplyLagrangianConstraints(double &chi_2, double &delta);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint)=0;
  virtual int numberOfEquations()=0;
  virtual bool isConverged();

 protected:
  virtual TVectorD Value(TVectorD &v)=0;
  virtual double ChiSquare(TMatrixT<double> delta_alpha,TMatrixT<double> lambda,TMatrixT<double> D,TMatrixT<double> d);
  virtual double ChiSquareUsingInitalPoint(TMatrixT<double> alpha,TMatrixT<double> lambda);
  virtual double ConstraintDelta(TVectorT<double> par);
  virtual TMatrixT<double> ComputeVariance();

  TVectorD par_first;
  TMatrixTSym<double> cov_first;
  TVectorD par;
  TVectorD par_prev;
  TMatrixTSym<double> cov;
  
  double chi2, chi2prev, delta;
  TVectorD maxStep;
  TVectorD w; // weight for constraints
  const MagneticField* field;
  bool debug;

 private:
  TMatrixT<double> Derivative();
  TVectorT<double> convertToVector(TMatrixT<double> M);
  TMatrixT<double> convertToMatrix(TVectorT<double> V);
  TVectorD WeightedValue(TVectorD &v){return weight_*Value(v);}
  double epsilon_, weight_,MaxDelta_;

  TMatrixTSym<double> V_alpha0_inv;
  TMatrixT<double> D;
  TMatrixTSym<double> V_D;
  double ScaleFactor;

  TMatrixT<double> V_corr_prev;
  
};
#endif
