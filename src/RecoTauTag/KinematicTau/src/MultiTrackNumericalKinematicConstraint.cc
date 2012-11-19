#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"

MultiTrackNumericalKinematicConstraint::MultiTrackNumericalKinematicConstraint(double weight):
  epsilon_(0.00001),
  weight_(weight)
{

}


bool MultiTrackNumericalKinematicConstraint::ApplyLagrangianConstraints(double &chi_2, double &delta){
  // get values
  AlgebraicVector val(numberOfEquations(),0);
  val=Value(par);  
  // obtain the derivatives
  AlgebraicMatrix g(numberOfEquations(),par.num_row(),0);
  g=Derivative();
  //copy covariance to sym. matrix
  AlgebraicSymMatrix in_cov_sym(cov_first.num_row(),0);
  for(int i = 1; i<=cov_first.num_row(); i++){
    for(int j = 1; j<=cov_first.num_row(); j++){
      if(i<=j) in_cov_sym(i,j) = cov_first(i,j);
    }
  }
  // delta_alpha
  AlgebraicVector delta_alpha=par_first-par;

  //debug code
  AlgebraicSymMatrix v_g_sym = in_cov_sym.similarity(g);
  int ifl1 = 0;
  v_g_sym.invert(ifl1);
  if(ifl1 !=0) {
    std::cout << "Fit failed: unable to invert SYM gain matrix\n" << std::endl;
    return false;
  }
  // v_g_sym now valid
  //full math case now!
  AlgebraicVector lambda = v_g_sym *(g*delta_alpha + val);
  //final parameters
  AlgebraicVector finPar = par -  in_cov_sym * g.T() * lambda;
  //covariance matrix business:
  AlgebraicMatrix mFactor = in_cov_sym *(v_g_sym.similarityT(g))* in_cov_sym;
  //refitted covariance
  AlgebraicMatrix rCov = in_cov_sym - mFactor;
  // copy new covariance
  // Force to be symmetric covariance (verify later)
  AlgebraicSymMatrix r_cov_sym(cov_first.num_row(),0);
  for(int i = 1; i<=cov_first.num_row(); i++){
    for(int j = 1; j<=cov_first.num_row(); j++){
      cov(i,j)=(rCov(i,j)+rCov(j,i))/2;
    }
  }
  //copy new par
  delta=0;
  for(int i = 1; i<=par.num_row(); i++){delta+=sqrt((par(i)-finPar(i))*(par(i)-finPar(i)));}
  par=finPar;

  AlgebraicVector chi  = lambda.T()*(g*delta_alpha  + val);
  chi_2=0;
  for(int i = 1; i<=chi.num_row(); i++){chi_2+=chi(i);}
  chi2=chi;  
  return true;
}


AlgebraicMatrix MultiTrackNumericalKinematicConstraint::Derivative(){
  AlgebraicMatrix d(numberOfEquations(),par.num_row(),0);
  AlgebraicVector par_plus(par.num_row(),0);
  AlgebraicVector value(numberOfEquations(),0);
  AlgebraicVector value_plus(numberOfEquations(),0);
  for(int j=1;j<=par.num_row();j++){
    for(int i=1;i<=par.num_row();i++){
      par_plus(i)=par(i);
      if(i==j) par_plus(i)=par(i)+epsilon_;
    }
    value=Value(par);
    value_plus=Value(par_plus);
    for(int i=1; i<=numberOfEquations();i++){
      d(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  return d;
}


