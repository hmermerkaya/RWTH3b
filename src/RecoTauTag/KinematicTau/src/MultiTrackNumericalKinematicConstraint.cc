#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TDecompBK.h"


MultiTrackNumericalKinematicConstraint::MultiTrackNumericalKinematicConstraint(bool debugflag,double weight):
  debug(debugflag),
  epsilon_(0.00001),
  weight_(weight)
{
  debugflag=false;
}


bool MultiTrackNumericalKinematicConstraint::ApplyLagrangianConstraints(double &chi_2, double &delta){
  //debug=false;
  //if(debug) std::cout << "MultiTrackNumericalKinematicConstraint::ApplyLagrangianConstraints start" <<std::endl;
  // Setup intial values
  TMatrixT<double> alpha_0=convertToMatrix(par);
  TMatrixT<double> alpha_A=convertToMatrix(par_first);
  TMatrixT<double> delta_alpha=alpha_0-alpha_A;
  TMatrixT<double> D=Derivative();
  TMatrixT<double> d=convertToMatrix(Value(par));

  //debug code
  TMatrixTSym<double> V_alpha0=cov;
  TMatrixTSym<double> V_D_inv=V_alpha0;
  V_D_inv.Similarity(D);
  double det = V_D_inv.Determinant();
  /*
  if(debug)std::cout << "===================================" << std::endl;
  if(debug)std::cout << "det " << det << std::endl;
  if(debug)std::cout << "V_alpha0" << std::endl;
  if(debug)V_alpha0.Print();
  if(debug)std::cout << "D" << std::endl;
  if(debug)D.Print();
  if(debug)std::cout << "VD_inverse" << std::endl;
  if(debug)V_D_inv.Print();
  if(debug)std::cout << "===================================" << std::endl;
  */
  if(fabs(det)<=1E-33){
    std::cout << "Fit failed: unable to invert SYM gain matrix " << det << " \n" << std::endl;
    return false;
  }
  TDecompBK Inverter(V_D_inv);
  TMatrixTSym<double> V_D=Inverter.Invert();

  // solve equations
  TMatrixT<double> lambda=V_D*(D*delta_alpha+d);
  TMatrixT<double> DT=D; DT.T();
  TMatrixT<double> alpha=alpha_0-V_alpha0*(DT*lambda);
  TVectorD finPar=convertToVector(alpha);

  // do while loop to see if the convergance criteria are satisfied
  double chi2prev(chi2), sf(1);//scale fraction
  bool MaxStepOk(true);
  TVectorT<double> dpar(par.GetNrows());

  for(int iter=0;iter<10;iter++){
    for(int i=0;i<par.GetNrows();i++) dpar(i)=sf*(finPar(i)-par(i)); 
    ///////////////////////////////////////////////////////////////////////////  
    // Check step size 
    MaxStepOk=true;
    double ratio(1.0);
    for(int i=0;i<par.GetNrows();i++){
      if(fabs(dpar(i))>maxStep(i)){ MaxStepOk=false; if(fabs(dpar(i)/maxStep(i))>ratio){ratio=1.1*fabs(dpar(i)/maxStep(i));} }
    }
    if(ratio>1)sf*=1.0/ratio;
    //std::cout << "par scale " << sf << " i= " << iter << std::endl;
    if(!MaxStepOk) continue;
    //////////////////////////////////////////////////////////////////////////
    // Compute convergence criteria.... chi2 and improvement in the constraints  
    // check that new chi square is an improvment
    TMatrixT<double> delta_alphaSf=delta_alpha+sf*(alpha-alpha_0);
    TMatrixT<double> lambdaT=lambda; lambdaT.T();
    TMatrixT<double> chisquare=lambdaT*(D*delta_alphaSf+d);
    chi2=0;
    for(int i = 0; i<chisquare.GetNrows(); i++){chi2+=chisquare(i,0);}
    chi_2=chi2;
    // check that the new constraint vector is an improvment
    TVectorD dfin(numberOfEquations());
    finPar=par+sf*dpar;
    dfin=Value(finPar);
    double d_C(0),dfin_C(0);
    for(int i = 0; i<d.GetNrows(); i++){
      d_C+=d(i,0)*d(i,0)*w(i);
      dfin_C+=dfin(i)*dfin(i)*w(i);
      //std::cout << "constriants " << dfin_C << " " << dfin(i) << " " << w(i) << std::endl;
    }
    std::cout << "constriants " << sf  << std::endl;
    delta=sqrt(dfin_C);
    dfin.Print();
    // correct step size if needed
    if(/*chi2<=chi2prev &&*/ delta<=sqrt(d_C)) break;
    sf*=0.5; 
    std::cout << "chi or conv scale " << sf << " i= " << iter << " dchi2 " << chi2-chi2prev << " deltla " << delta-sqrt(d_C) <<  std::endl;
  }
  //correct finPar to new stepsize
  finPar=par+sf*dpar;
  //covariance matrix with modified step size
  TMatrixTSym<double> DTV_DD=V_D.SimilarityT(D);
  TMatrixT<double> DTV_DDV=DTV_DD*V_alpha0;
  TMatrixT<double> VDTV_DDV=V_alpha0*DTV_DDV;
  TMatrixT<double> CovCor=VDTV_DDV;
  TMatrixT<double> V_alpha = V_alpha0-CovCor;
  //copy new par
  std::cout << "alpha_0" << std::endl;  
  par.Print();
  par=finPar;
  for(int i=0; i<cov.GetNrows();i++){
    for(int j=0; j<=i;j++){
      cov(i,j)=V_alpha(i,j);
    }
  }
  std::cout << "alpha" << std::endl;
  par.Print();
  std::cout << " delta " << delta << std::endl; 
  //if(debug)std::cout << "MultiTrackNumericalKinematicConstraint::ApplyLagrangianConstraints end" <<std::endl;
  return true;
}


TMatrixD MultiTrackNumericalKinematicConstraint::Derivative(){
  //std::cout << "MultiTrackNumericalKinematicConstraint::Derivative start" << std::endl;
  bool hasdebug=debug;
  debug=false;
  TMatrixD D(numberOfEquations(),par.GetNrows());
  TVectorD par_plus(par.GetNrows());
  TVectorD value(numberOfEquations());
  TVectorD value_plus(numberOfEquations());
  for(int j=0;j<par.GetNrows();j++){
    for(int i=0;i<par.GetNrows();i++){
      par_plus(i)=par(i);
      if(i==j) par_plus(i)=par(i)+epsilon_;
    }
    value=Value(par);
    value_plus=Value(par_plus);
    for(int i=0; i<numberOfEquations();i++){
      D(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  debug=hasdebug;
  //std::cout << "MultiTrackNumericalKinematicConstraint::Derivative end" << std::endl;
  return D;
}


bool MultiTrackNumericalKinematicConstraint::isConverged(double chi_2, double delta){
  if(delta<0.01) return true;
  return false;
}



TVectorT<double> MultiTrackNumericalKinematicConstraint::convertToVector(TMatrixT<double> M){
  TVectorT<double> V(M.GetNrows()); 
  for(int i=0; i<M.GetNrows();i++){
    V(i)=M(i,0);
  }
  return V;
}


TMatrixT<double> MultiTrackNumericalKinematicConstraint::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i<M.GetNrows();i++){
    M(i,0)=V(i);
  }
  return M;
}
