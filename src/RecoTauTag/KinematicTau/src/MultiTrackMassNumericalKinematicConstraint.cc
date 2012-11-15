#include "RecoTauTag/KinematicTau/interface/MultiTrackMassNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"

AlgebraicVector MultiTrackMassNumericalKinematicConstraint::Compute(const std::vector<KinematicState> states,const GlobalPoint& point, const std::vector<double> &epsilonPar, const std::vector<double> &epsilonPos) const{
  if(states.size()<nPart) throw VertexException("MultiTrackMassNumericalKinematicConstraint::not enough states given");
  double sumEnergy(0),sumPx(0),sumPy(0),sumPz(0),m(0),px(0),py(0),pz(0),a(0);
  for (unsigned int i=0;i<nPart;++i) {
    a = -states.at(i).particleCharge() * states.at(i).magneticField()->inInverseGeV(states.at(i).globalPosition()).z();
    
    px = states.at(i).kinematicParameters()(3)+epsilonPar.at(3+i*npardim);
    py = states.at(i).kinematicParameters()(4)+epsilonPar.at(4+i*npardim);
    pz = states.at(i).kinematicParameters()(5)+epsilonPar.at(5+i*npardim);
    m  = states.at(i).kinematicParameters()(6)+epsilonPar.at(6+i*npardim);
   
    sumEnergy  += sqrt(m*m+px*px+py*py+pz*pz);
    sumPx      += px- a*((point.y()+epsilonPos.at(1)) - (states.at(i).kinematicParameters()(1)+epsilonPar.at(1+i*npardim)));
    sumPy      += py+ a*((point.x()+epsilonPos.at(0)) - (states.at(i).kinematicParameters()(0)+epsilonPar.at(0+i*npardim)));
    sumPz      += pz;
  }
  double j_m = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz;
  AlgebraicVector res(1,0);
  res(1)  = (sumEnergy*sumEnergy - j_m - mass*mass)*GetWeight();
  std::cout << "Mass " << res(1) << std::endl;
  return res;
}
