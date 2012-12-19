#ifndef ThreeProngTauSolver_h
#define ThreeProngTauSolver_h

#include "TVector3.h" 
#include "TLorentzVector.h" 

class  ThreeProngTauSolver {
 public:
  // constructor and Destructor
  ThreeProngTauSolver(){};
  ~ThreeProngTauSolver(){};

  void quadratic(double &x_plus,double &x_minus,double a, double b, double c);
  void AnalyticESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1);
  void NumericalESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1);
  void SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2,bool rotateback=true);
  bool SetTauDirectionatThetaGJMax(TLorentzVector a1, double &theta,double &phi);  
};

#endif
