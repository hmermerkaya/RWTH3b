#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"

MultiTrackNumericalKinematicConstraint::MultiTrackNumericalKinematicConstraint(double weight):
  npardim(7),
  nposdim(3),
  distance_epsilon(0.000001),
  momentum_epsilon(0.00001),
  weight_(weight)
{
}

AlgebraicVector  MultiTrackNumericalKinematicConstraint::value(const std::vector<KinematicState> states,
                        const GlobalPoint& point) const
{
 if(states.size()<2) throw VertexException("MultiTrackNumericalKinematicConstraint::value not enough states given");
 AlgebraicVector      value(numberOfEquations(),0);
 value=ComputeDefault(states,point);
 return  ComputeDefault(states,point);
}

AlgebraicMatrix MultiTrackNumericalKinematicConstraint::parametersDerivative(const std::vector<KinematicState> states,
                                      const GlobalPoint& point) const
{
  if(states.size()<2) throw VertexException("MultiTrackNumericalKinematicConstraint::parametersDerivative not enough states given");
  AlgebraicMatrix      res(numberOfEquations(),states.size()*npardim,0);
  AlgebraicVector      value(numberOfEquations(),0);
  AlgebraicVector      valueplus(numberOfEquations(),0);
  std::vector<double>  epsilonPar(npardim*states.size(),0);
  std::vector<double>  epsilonPos(nposdim,0);
  for (unsigned int j=0;j<states.size();j++) {
    for(unsigned int k=1;k<=npardim;k++){
      value=ComputeDefault(states,point);
      int idx=k+j*npardim-1;
      if(k<=3) epsilonPar.at(idx)=distance_epsilon;
      else epsilonPar.at(idx)=momentum_epsilon;
      if(epsilonPar.at(idx)!=0){
	valueplus=Compute(states,point,epsilonPar,epsilonPos);
	for(unsigned int i=1;i<=(unsigned int)numberOfEquations();i++){
	  res(i,k+j*npardim) = (valueplus(i)-value(i))/epsilonPar.at(idx); 
	}
      }
      else{
	VertexException("MultiTrackNumericalKinematicConstraint::parametersDerivative Epsilon=0 ");
      }
      epsilonPar.at(idx)=0;
    }
  }
  return res;
}

AlgebraicMatrix MultiTrackNumericalKinematicConstraint::positionDerivative(const std::vector<KinematicState> states,
								      const GlobalPoint& point) const
{
  if(states.size()<2) throw VertexException("MultiTrackNumericalKinematicConstraint::positionDerivative not enough states given");
  AlgebraicMatrix      res(numberOfEquations(),nposdim,0);
  AlgebraicVector      value(numberOfEquations(),0);
  AlgebraicVector      valueplus(numberOfEquations(),0);
  std::vector<double>  epsilonPar(npardim*states.size(),0);
  std::vector<double>  epsilonPos(nposdim,0);

  for(unsigned int k=1;k<=nposdim;k++){
    value=ComputeDefault(states,point);
    int idx=k-1;
    epsilonPos.at(idx)=distance_epsilon;
    if(epsilonPos.at(idx)!=0){
      valueplus=Compute(states,point,epsilonPar,epsilonPos);
      for(unsigned int i=1;i<=(unsigned int)numberOfEquations();i++){
	res(i,k) = (valueplus(i)-value(i))/epsilonPos.at(idx);
      }
    }
    else{
      VertexException("MultiTrackNumericalKinematicConstraint::positionDerivative Epsilon=0 ");
    }
    epsilonPos.at(idx)=0;
  }
  return res;
}

AlgebraicVector  MultiTrackNumericalKinematicConstraint::ComputeDefault(const std::vector<KinematicState> states,const GlobalPoint& point) const{
  std::vector<double>  epsilonPar(npardim*states.size(),0);
  std::vector<double>  epsilonPos(nposdim,0);
  return Compute(states,point,epsilonPar,epsilonPos);
}
