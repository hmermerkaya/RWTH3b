#include "RecoTauTag/KinematicTau/interface/MultiTrackSmartPointingNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"

AlgebraicVector MultiTrackSmartPointingNumericalKinematicConstraint::Compute(const std::vector<KinematicState> states,const GlobalPoint& point, const std::vector<double> &epsilonPar, const std::vector<double> &epsilonPos) const{
  if(states.size()<2) throw VertexException("MultiTrackSmartPointingNumericalKinematicConstraint::not enough states given");
  int num = states.size();
  if(num<2) throw VertexException("MultiTrackPointingKinematicConstraint::value <2 states passed");

  //2 equations (for all tracks)                                                                                                                                                                                                       
  AlgebraicVector  vl(2,0);
  double dx = (point.x()+epsilonPos.at(0)) - refPoint.x();
  double dy = (point.y()+epsilonPos.at(1)) - refPoint.y();
  double dz = (point.z()+epsilonPos.at(2)) - refPoint.z();
  
  double sumPx(0),sumPy(0),sumPz(0),px(0),py(0),pz(0),a(0);
  for(unsigned int i=0;i<states.size();i++){
    //std::cout << "state " << i << " charge " << states.at(i).particleCharge() << " m " << states.at(i).kinematicParameters()(6) <<  std::endl;
    a  = -states.at(i).particleCharge() * states.at(i).magneticField()->inInverseGeV(states.at(i).globalPosition()).z();
    px = states.at(i).kinematicParameters()(3)+epsilonPar.at(3+i*npardim);
    py = states.at(i).kinematicParameters()(4)+epsilonPar.at(4+i*npardim);
    pz = states.at(i).kinematicParameters()(5)+epsilonPar.at(5+i*npardim);
    
    sumPx      += px;// - a*((point.y()+epsilonPos.at(1)) - (states.at(i).kinematicParameters()(1)+epsilonPar.at(1+i*npardim)));
    sumPy      += py;// + a*((point.x()+epsilonPos.at(0)) - (states.at(i).kinematicParameters()(0)+epsilonPar.at(0+i*npardim)));
    sumPz      += pz;
  }
  
  /*
  double dT = sqrt(pow(dx,2) + pow(dy,2));
  double ds = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
  double pT = sqrt(pow(sumPx,2) + pow(sumPy,2));
  double pSum = sqrt(pow(sumPx,2) + pow(sumPy,2) + pow(sumPz,2));
  
  vl(1) = (dT - dx)/dy + (sumPx - pT)/sumPy;
  vl(2) = (ds - dT)/dz + (pT - pSum)/sumPz;


  vl(1)*=1000;
  vl(2)*=1000;
  */
  
  double DTdotPT = sumPx*dx+sumPy*dy;
  double DTdotDT = dx*dx+dy*dy;
  double PTdotPT = sumPx*sumPx+sumPy*sumPy;

  double DdotP = DTdotPT+sumPz*dz;
  double DdotD = DTdotDT+dz*dz;
  double PdotP = PTdotPT+sumPz*sumPz;

  double cosphi2=(DTdotPT*DTdotPT)/(DTdotDT*PTdotPT);
  double costheta2=(DdotP*DdotP)/(DdotD*PdotP);

  vl(1) = GetWeight()*(1-cosphi2);
  if(DTdotPT<0) vl(1) = GetWeight()*(1+cosphi2);
  vl(2) = GetWeight()*(1-costheta2);
  if(DdotP<0) vl(2) = GetWeight()*(1+costheta2);
  //  std::cout << "dx: " << dx/sqrt(DdotD) << " dy " << dy/sqrt(DdotD) << " dz " << dz/sqrt(DdotD) << " px " << sumPx/sqrt(PdotP) << " py " << sumPy/sqrt(PdotP) << " pz " << sumPz/sqrt(PdotP) << " " << cosphi2 << " " << costheta2 << std::endl;   
  //std::cout << "v1: " << vl(1) << " v2: " << vl(2) << " " << GetWeight() << std::endl;
  
  return vl;
}


