#include "RecoTauTag/KinematicTau/interface/TauA1NuNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertexFactory.h"
#include "TLorentzVector.h"

AlgebraicVector TauA1NuNumericalKinematicConstraint::Value(AlgebraicVector &v){
  TLorentzVector a1(v(4),v(5),v(6), v(7));
  TLorentzVector nu(v(8),v(9),v(10),sqrt(v(8)*v(8)+v(9)*v(9)+v(10)*v(10)));

  TVector3 sv(v(1),v(2),v(3));
  TVector3 pv=pv_inital;
  TVector3 TauDir=sv-pv;
  double phi(TauDir.Phi()),theta(TauDir.Theta());
  TLorentzVector a1rot=a1;
  a1rot.RotateZ(-phi);
  a1rot.RotateY(-theta);
  TLorentzVector nurot=nu;
  nurot.RotateZ(-phi);
  nurot.RotateY(-theta);

  TLorentzVector nurotfixed(a1rot.Px(),a1rot.Py(),nurot.Pz(),sqrt(a1rot.Pt()*a1rot.Pt()+nurot.Pz()*nurot.Pz()));
  AlgebraicVector res(3,0);
  TLorentzVector tau=a1rot+nurotfixed;
  res(1) = tau.M(); // mass constraint fixed to only float Pz for nu (ie only one value with huge errors per constraint) 
  res(2) = a1rot.Pt()-nurot.Pt(); // |Pt| balance constraint
  res(3) = 1+(a1rot.Px()*nurot.Px()+a1rot.Py()*nurot.Py())/(a1rot.Pt()*nurot.Pt()); // phi' constraint (back-to-back) 
  return res;
}


bool TauA1NuNumericalKinematicConstraint::ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint){
  AlgebraicVector inpar(10,0);
  AlgebraicMatrix incov(10,10,0);
  bool hasa1(false),hasnu(false);
  if(inStates.size()!=2) return false;
  for(unsigned int i=0; i<inStates.size();i++){
    AlgebraicVector    inPar = asHepVector<7>(inStates.at(i).kinematicParameters().vector());
    AlgebraicSymMatrix inCov = asHepMatrix<7>(inStates.at(i).kinematicParametersError().matrix());
    if(inStates.at(i).particleCharge()!=0){
      hasa1=true;
      // Parameter
      for(int j = 1; j<8; j++){inpar(j) = inPar(j);}
      
      // Error Matrix
      for(int j = 1; j<8; j++){
	for(int k = 1; k<8; k++){
	  incov(j,k)=inCov(j,k);
	}
      }
    }
    else{
      hasnu=true;
      // Parameter           
      for(int j = 4; j<7; j++){inpar(j+4) = inPar(j);}
      
      // Error Matrix
      for(int j = 4; j<7; j++){
        for(int k = 4; k<7; k++){
          incov(j+4,k+4)=inCov(j,k);
        }
      }
    }
  }
  if(hasa1 && hasnu){
    par=inpar;
    cov=incov;
    par_first=par;
    cov_first=cov;
    field=inStates.front().magneticField();
    return true;
  }
  return false;
}

std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > TauA1NuNumericalKinematicConstraint::ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint){
  std::cout << "r" << inStates.size() << std::endl;
  AlgebraicSymMatrix r_cov_sym(cov_first.num_row(),0);
  for(int i = 1; i<=cov_first.num_row(); i++){
    for(int j = 1; j<=cov_first.num_row(); j++){
      if(i<=j){
	r_cov_sym(i,j)  = cov(i,j);
      }
    }
  }
  AlgebraicSymMatrix pCov = r_cov_sym.sub(1,3);
  std::cout << "s" << inStates.size() << std::endl;
  //making resulting vertex
  float ndf=3.0;
  GlobalPoint vPos (par(1),par(2),par(3));
  VertexState st(vPos,GlobalError( asSMatrix<3>(pCov)));
  KinematicVertexFactory vFactory;
  RefCountedKinematicVertex rVtx = vFactory.vertex(st,chi2(1),ndf);

  //making refitted states of Kinematic Particles
  std::vector<KinematicState> ns;

  AlgebraicVector7 newParA1;
  AlgebraicVector7 newParNu;
  for(int i =1; i<8; i++){newParA1(i) = par(i);}
  for(int i =1; i<3; i++){newParNu(i) = par(i);}
  for(int i =8; i<11; i++){newParNu(i-4) = par(i);}
  // newParNu(7)=0;
  std::cout << "t" << inStates.size() << std::endl;

  AlgebraicMatrix ns_cov(14+3,14+3,0);
  for(unsigned int i=0; i<inStates.size();i++){
    if(inStates.at(i).particleCharge()!=0){
      std::cout <<"a" << std::endl;
      AlgebraicSymMatrix nCovariance = r_cov_sym.sub(1,7);
      TrackCharge chl = inStates.at(0).particleCharge();
      KinematicParameters nrPar(newParA1);
      KinematicParametersError nrEr(asSMatrix<7>(nCovariance));
      KinematicState newState(nrPar,nrEr,chl, field);
      ns.push_back(newState);
      std::cout <<"b" << std::endl;
      for(int i = 1; i<=cov.num_row(); i++){
	for(int j = 1; j<=cov.num_row(); j++){
	  ns_cov(i,j)=cov(i,j);
	}
      }
      std::cout <<"c" << std::endl;
    }
    else{
      AlgebraicSymMatrix nCovariance = r_cov_sym.sub(1,7);
      /*      for(int i = 4; i<=cov_first.num_row(); i++){
	for(int j = 4; j<=cov_first.num_row(); j++){
	  if(i<=j){nCovariance(i,j)=cov(i+4,j+4);}
	  if(j==7 ){
	    if(i!=j){nCovariance(i,j)=0;}
	    else{nCovariance(i,j)=0.001;}
	  }
	}
	}*/
      std::cout <<"d" << std::endl;
      for(int i = 1; i<=cov.num_row(); i++){
        for(int j = 1; j<=cov.num_row(); j++){
          ns_cov(j+7,i+7)=cov(i,j);
        }
      }
      for(int i = 1; i<=3; i++){
        for(int j = 1; j<=3; j++){
	  ns_cov(j+14,i+14)=cov(i,j);
	}
      }
      std::cout <<"e" << std::endl;
      TrackCharge chl = inStates.at(1).particleCharge();
      std::cout <<"f" << std::endl;
      KinematicParameters nrPar(newParNu);
      std::cout <<"g" << std::endl;
      KinematicParametersError nrEr(asSMatrix<7>(nCovariance));
      std::cout <<"h" << std::endl;
      KinematicState newState(nrPar,nrEr,chl,field);
      std::cout <<"i" << std::endl;
      ns.push_back(newState);
      std::cout <<"j" << std::endl;
    }
  }
  std::cout << "v" << inStates.size() << std::endl;
  std::cout << ns_cov << std::endl; 

  std::pair<std::vector<KinematicState>, AlgebraicMatrix> ns_m(ns,ns_cov);
  return std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix>, RefCountedKinematicVertex >(ns_m,rVtx);
}
