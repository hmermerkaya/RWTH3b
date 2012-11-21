#include "RecoTauTag/KinematicTau/interface/TauA1NuNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertexFactory.h"
#include "TLorentzVector.h"

AlgebraicVector TauA1NuNumericalKinematicConstraint::Value(AlgebraicVector &v){
  std::cout << v << std::endl;

  TLorentzVector a1(v(4),v(5),v(6),sqrt(v(4)*v(4)+v(5)*v(5)+v(6)*v(6)+v(7)*v(7)));
  TLorentzVector nu(v(8),v(9),v(10),sqrt(v(8)*v(8)+v(9)*v(9)+v(10)*v(10)));
  TLorentzVector a1_d=a1;
  std::cout<<"nu Phi  " <<nu.Phi()<< " nu Px " <<  nu.Px() << " nu Px " << nu.Py() << " nu Pt " <<  nu.Pz() <<std::endl;
  std::cout<<"a1 Phi  " <<a1.Phi()<< " a1 Px " <<  a1.Px() << " a1 Px " << a1.Py() << " a1 Pt " <<  a1.Pz() << " mass " << a1.M() << " E " << a1.E() <<std::endl;


  TVector3 sv(v(1),v(2),v(3));
  TVector3 pv=pv_inital;
  TVector3 TauDir=sv-pv;
  double phi(TauDir.Phi()),theta(TauDir.Theta());
  a1.RotateZ(-phi);
  a1.RotateY(-theta);
  nu.RotateZ(-phi);
  nu.RotateY(-theta);

  TLorentzVector nufixed(-a1.Px(),-a1.Py(),nu.Pz(),sqrt(a1.Pt()*a1.Pt()+nu.Pz()*nu.Pz()));
  AlgebraicVector res(3,0);
  TLorentzVector tau=a1+nufixed;
  res(1) = tau.M()-mtau_c; // mass constraint fixed to only float Pz for nu (ie only one value with huge errors per constraint) 
  res(2) = a1.Pt()-nu.Pt(); // |Pt| balance constraint
  res(3) = 1+(a1.Px()*nu.Px()+a1.Py()*nu.Py())/(a1.Pt()*nu.Pt()); // phi' constraint (back-to-back) 

  if(debug){
    std::cout << "------------>"<<std::endl;
    std::cout << "nufixed Phi  " <<nufixed.Phi()<< " nufixed Px " <<  nufixed.Px() << " nufixed Px " << nufixed.Py() << " nufixed Pt " <<  nufixed.Pt() <<std::endl;
    std::cout << "a1 Phi  " <<a1.Phi()<< " a1 Px " <<  a1.Px() << " a1 Px " << a1.Py() << " a1 Pt " <<  a1.Pt() << " mass " << a1.M()  <<std::endl;
    std::cout << "Tau  M" << tau.M() << " E " << tau.E() <<  " Px" << tau.Px() << " Py" << tau.Py() << " Pz " << tau.Pz() << std::endl;
    if(GenPart.isValid()){
      double TauMatchingDR_=0.4;
      for(reco::GenParticleCollection::const_iterator itr = GenPart->begin(); itr!= GenPart->end(); ++itr){
	if(itr->pdgId()==15){
	  const reco::GenParticle mytau=(*itr);
	  TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	  if(a1_d.DeltaR(mc)<TauMatchingDR_){
	    for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){
	      const reco::Candidate *dau=(itr)->daughter(i);
	      if(fabs(dau->pdgId())==20213){
		TVector3 TruthPvtx((itr)->vx(),(itr)->vy(),(itr)->vz());
		TVector3 TruthSvtx(dau->vx(),dau->vy(),dau->vz());
		TVector3 TruthDir=TruthSvtx-TruthPvtx; 
		std::cout << " TauDir Phi " << TauDir.Phi() << " Theta " << TauDir.Theta() << " Mag " << TauDir.Mag() << " E " << tau.E() << std::endl;
		std::cout << " Truth Dir Phi " << TruthDir.Phi() << " Theta " << TruthDir.Theta() << " Mag " << TruthDir.Mag() << " E " << mc.E() << std::endl;
		std::cout << "PVT Truth Vx"  << TruthPvtx.Px() << " Vy " << TruthPvtx.Py() <<  " Vz " << TruthPvtx.Pz() << " " << std::endl;
		std::cout << "PV Current Vx"  << pv.Px() << " Vy " << pv.Py() <<  " Vz " << pv.Pz() << " " << std::endl;
		std::cout << "SVT Truth Vx"  << TruthSvtx.Px() << " Vy " << TruthSvtx.Py() <<  " Vz " << TruthSvtx.Pz() << " " << std::endl;
		std::cout << "SV Current Vx"  << sv.Px() << " Vy " << sv.Py() <<  " Vz " << sv.Pz() << " " << std::endl;
		
	      }
	    }
	  }
	}
      }
    }
  }
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
  AlgebraicSymMatrix r_cov_sym(cov_first.num_row(),0);
  for(int i = 1; i<=cov_first.num_row(); i++){
    for(int j = 1; j<=cov_first.num_row(); j++){
      if(i<=j){
	r_cov_sym(i,j)  = cov(i,j);
      }
    }
  }
  AlgebraicSymMatrix pCov = r_cov_sym.sub(1,3);
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

  AlgebraicMatrix ns_cov(14+3,14+3,0);
  for(unsigned int i=0; i<inStates.size();i++){
    if(inStates.at(i).particleCharge()!=0){
      AlgebraicSymMatrix nCovariance = r_cov_sym.sub(1,7);
      TrackCharge chl = inStates.at(0).particleCharge();
      KinematicParameters nrPar(newParA1);
      KinematicParametersError nrEr(asSMatrix<7>(nCovariance));
      KinematicState newState(nrPar,nrEr,chl, field);
      ns.push_back(newState);
      for(int i = 1; i<=cov.num_row(); i++){
	for(int j = 1; j<=cov.num_row(); j++){
	  ns_cov(i,j)=cov(i,j);
	}
      }
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
      TrackCharge chl = inStates.at(1).particleCharge();
      KinematicParameters nrPar(newParNu);
      KinematicParametersError nrEr(asSMatrix<7>(nCovariance));
      KinematicState newState(nrPar,nrEr,chl,field);
      ns.push_back(newState);
    }
  }

  std::pair<std::vector<KinematicState>, AlgebraicMatrix> ns_m(ns,ns_cov);
  return std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix>, RefCountedKinematicVertex >(ns_m,rVtx);
}
