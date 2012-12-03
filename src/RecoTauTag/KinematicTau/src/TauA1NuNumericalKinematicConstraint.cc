#include "RecoTauTag/KinematicTau/interface/TauA1NuNumericalKinematicConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertexFactory.h"
#include "TLorentzVector.h"

TauA1NuNumericalKinematicConstraint::TauA1NuNumericalKinematicConstraint(const reco::Vertex &primaryVertex,double mtau,edm::Handle<reco::GenParticleCollection> &GenPart_,double weight,bool debug_):
  MultiTrackNumericalKinematicConstraint(debug_),
  pv_inital(primaryVertex),
  mtau_c(mtau),
  GenPart(GenPart_)
{
  pv_inital=primaryVertex;
  if(!GenPart.isValid())debug=false;
  TVectorD maxStepTemp(npar);
  maxStepTemp(tau_phi)=0.001;
  maxStepTemp(tau_theta)=0.001;
  maxStepTemp(a1_px)=0.5;
  maxStepTemp(a1_py)=0.5;
  maxStepTemp(a1_pz)=0.5;
  maxStepTemp(a1_m)=0.1;
  maxStepTemp(nu_px)=2.0;
  maxStepTemp(nu_py)=2.0;
  maxStepTemp(nu_pz)=2.0;
  maxStep.ResizeTo(maxStepTemp);
  maxStep=maxStepTemp;
}



TVectorD TauA1NuNumericalKinematicConstraint::Value(TVectorD &v){
  TLorentzVector a1(v(a1_px),v(a1_py),v(a1_pz),sqrt(v(a1_m)*v(a1_m)+v(a1_px)*v(a1_px)+v(a1_py)*v(a1_py)+v(a1_pz)*v(a1_pz)));
  TLorentzVector nu(v(nu_px),v(nu_py),v(nu_pz),sqrt(v(nu_px)*v(nu_px)+v(nu_py)*v(nu_py)+v(nu_pz)*v(nu_pz)));
  TLorentzVector a1_d=a1;
  TLorentzVector nu_d=nu;
  double phi(v(tau_phi)),theta(v(tau_theta));
  a1.RotateZ(-phi);
  a1.RotateY(-theta);
  nu.RotateZ(-phi);
  nu.RotateY(-theta);

  TLorentzVector nufixed=nu;//(-a1.Px(),-a1.Py(),nu.Pz(),sqrt(a1.Pt()*a1.Pt()+nu.Pz()*nu.Pz()));
  TVectorD res(3);
  TLorentzVector tau=a1+nufixed;
  res(0) = tau.M2()-mtau_c*mtau_c;                                // mass constraint fixed to only float Pz for nu (ie only one value with huge errors per constraint) 
  res(1) = a1.Px()+nu.Px(); //a1.Pt()-nu.Pt();                    // |Pt| balance constraint
  res(2) = a1.Py()+nu.Py(); //1+(a1.Px()*nu.Px()+a1.Py()*nu.Py()) // (a1.Pt()*nu.Pt()); // phi' constraint (back-to-back) 
  /*
  if(debug){
    std::cout << "------------>"<<std::endl;
    std::cout << "Parameters: " << std::endl; v.Print();
    std::cout << "Constraints: " << std::endl; res.Print();
    std::cout << "nufixed Phi  " << nufixed.Phi()<< " nufixed Px " <<  nufixed.Px() << " nufixed Px " << nufixed.Py() << " nufixed Pt " <<  nufixed.Pt() <<std::endl;
    std::cout << "a1 Phi  " <<a1.Phi()<< " a1 Px " <<  a1.Px() << " a1 Px " << a1.Py() << " a1 Pt " <<  a1.Pt() << " mass " << a1.M()  <<std::endl;
    std::cout << "Tau  M" << tau.M() << " E " << tau.E() <<  " Px" << tau.Px() << " Py" << tau.Py() << " Pz " << tau.Pz() << std::endl;
    std::cout << "------------>"<<std::endl;
    if(GenPart.isValid()){
      double TauMatchingDR_=0.4;
      for(reco::GenParticleCollection::const_iterator itr = GenPart->begin(); itr!= GenPart->end(); ++itr){
	if(itr->pdgId()==15){
	  const reco::GenParticle mytau=(*itr);
	  TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	  if(a1_d.DeltaR(mc)<TauMatchingDR_){
	    TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	    for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){
	      const reco::Candidate *dau=(itr)->daughter(i);
	      if(fabs(dau->pdgId())==20213){
		// Print info from matching tau
		std::cout << "Tau Truth Px " << mc.Px() << " Py " << mc.Py() << " Pz " << mc.Pz() << " E " << mc.E() << " M " << mc.M() << std::endl;
		tau.RotateY(theta);
		tau.RotateZ(phi);
		std::cout << "Tau Truth Px " << tau.Px() << " Py " << tau.Py() << " Pz " << tau.Pz() << " E " << tau.E() << " M " << tau.M() << std::endl;
		// print a1 info
		std::cout << "A1 Truth Px " << dau->p4().Px() << " Py " << dau->p4().Py() << " Pz " << dau->p4().Pz() << " E " << dau->p4().E() << " M " << dau->p4().M() <<  std::endl;
		a1.RotateY(theta);
                a1.RotateZ(phi);
		std::cout << "A1 d Px " << a1_d.Px() << " Py " << a1_d.Py() << " Pz " << a1_d.Pz() << " E " << a1_d.E() << " M " << a1_d.M() << std::endl;
		TVector3 TruthPvtx((itr)->vx(),(itr)->vy(),(itr)->vz());
		TVector3 TruthSvtx(dau->vx(),dau->vy(),dau->vz());
  */
  /*
		TVector3 TruthDir=TruthSvtx-TruthPvtx; 
		std::cout << " TauDir Phi "    << phi            << " Theta " << theta << " Mag " << TauDir.Mag() << " E " << tau.E() << std::endl;
		std::cout << " Truth Dir Phi " << TruthDir.Phi() << " Theta " << TruthDir.Theta() << " Mag " << TruthDir.Mag() << " E " << mc.E() << std::endl;
		std::cout << "PVT Truth Vx"    << TruthPvtx.Px() << " Vy " << TruthPvtx.Py() <<  " Vz " << TruthPvtx.Pz() << " " << std::endl;
		std::cout << "PV Current Vx"   << pv.Px()        << " Vy " << pv.Py() <<  " Vz " << pv.Pz() << " " << std::endl;
		std::cout << "SVT Truth Vx"    << TruthSvtx.Px() << " Vy " << TruthSvtx.Py() <<  " Vz " << TruthSvtx.Pz() << " " << std::endl;
		std::cout << "SV Current Vx"   << sv.Px()        << " Vy " << sv.Py() <<  " Vz " << sv.Pz() << " " << std::endl;
		*/
  /*
	      }
	      if(fabs(dau->pdgId())==16){
                nu.RotateY(theta);
		nu.RotateZ(phi);
		std::cout << "Neutrino Truth Px " << dau->p4().Px() << " Py " << dau->p4().Py() << " Pz " << dau->p4().Pz() << " E " << dau->p4().E() << std::endl;
		//std::cout << "Neutrino    Px " << nu.Px() << " Py " << nu.Py() << " Pz " << nu.Pz() << " E " << nu.E() << std::endl;
		std::cout << "Neutrino d  Px " << nu_d.Px() << " Py " << nu_d.Py() << " Pz " << nu_d.Pz() << " E " << nu_d.E() << std::endl;
	      }
	    }
	  }
	}
      }
    }
  }*/
  return res;
}


bool TauA1NuNumericalKinematicConstraint::ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint){
  //if(debug) std::cout << "TauA1NuNumericalKinematicConstraint::ConfigureIntialState start" << std::endl;
  // setup 13 by 13 matrix
  TVectorD inpar_sv(npar);
  TVectorD inpar(npar);
  TMatrixTSym<double> incov(npar);
  TMatrixTSym<double> incov_temp(norigpar);
  bool hasa1(false),hasnu(false);
  if(inStates.size()!=2) return false;
  for(unsigned int i=0; i<inStates.size();i++){
    AlgebraicVector    inPar = asHepVector<7>(inStates.at(i).kinematicParameters().vector());
    AlgebraicSymMatrix inCov = asHepMatrix<7>(inStates.at(i).kinematicParametersError().matrix());
    if(inStates.at(i).particleCharge()!=0){
      hasa1=true;
      // Parameter
      for(int j=0; j<npardim; j++){
	if(j>par_vz){inpar(j+(a1_px-par_px)) = inPar(j+1);}
	else{inpar_sv(j)=inPar(j+1);}
      }
      // Error Matrix
      for(int j=0; j<npardim; j++){
	for(int k=0; k<npardim; k++){
	  incov_temp(j+nposdim,k+nposdim)=inCov(j+1,k+1);
	}
      }
    }
    else{
      hasnu=true;
      // Parameter     
      for(int j = par_px; j<=par_pz; j++){inpar(j+(nu_px-par_px))=inPar(j+1);}
      // Error Matrix
      for(int j = par_px; j<=par_pz; j++){
        for(int k = par_px; k<=par_pz; k++){
          incov_temp(j+(nposdim+npardim-par_px),k+(nposdim+npardim-par_px))=inCov(j+1,k+1);
        }
      }
    }
  }
  // Get vertices and convert into parameters and cov
  pv=TVector3(pv_inital.x(),pv_inital.y(),pv_inital.z());
  sv=TVector3(inpar_sv(par_vx),inpar_sv(par_vy),inpar_sv(par_vz));
  TauDir=sv-pv;
  math::Error<nposdim>::type pvCov;
  pv_inital.fill(pvCov);
  for(int i = 0; i<nposdim; i++){
    for(int j = 0; j<=i; j++){
      incov_temp(i,j)=pvCov(i,j);
      incov_temp(j,i)=pvCov(i,j);
    }
  }

  inpar(tau_phi)=TauDir.Phi();
  inpar(tau_theta)=TauDir.Theta();

  // now transform errors using error propogation
  double xs(sv.X()),ys(sv.Y()),zs(sv.Z()),xp(pv.X()),yp(pv.Y()),zp(pv.Z());
  
  TMatrixT<double> D(norigpar,npar);
  D(0,0)=-(-pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.1e1 / 0.2e1) - (double) (xs - xp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.3e1 / 0.2e1) * (double) (-2 * xs + 2 * xp) / 0.2e1) * pow((double) (1 - (int) pow((double) (xs - xp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp)), -0.1e1 / 0.2e1); //dphi/dxp
  D(1,0)=(double) (xs - xp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.3e1 / 0.2e1) * (double) (2 * ys - 2 * yp) * pow((double) (1 - (int) pow((double) (xs - xp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp)), -0.1e1 / 0.2e1) / 0.2e1; //dphi/dyp 
  D(2,0)=0; //dphi/dzp 
  D(3,0)=-(pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.1e1 / 0.2e1) - (double) (xs - xp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.3e1 / 0.2e1) * (double) (2 * xs - 2 * xp) / 0.2e1) * pow((double) (1 - (int) pow((double) (xs - xp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp)), -0.1e1 / 0.2e1); //dphi/dxs 
  D(4,0)= (double) (xs - xp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp), -0.3e1 / 0.2e1) * (double) (2 * ys - 2 * yp) * pow((double) (1 - (int) pow((double) (xs - xp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp)), -0.1e1 / 0.2e1) / 0.2e1; //dphi/dys 
  D(5,0)=0; //dphi/dzs 

  D(0,1)=(double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (-2 * xs + 2 * xp) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1) / 0.2e1; //dtheta/dxp           
  D(1,1)= (double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (2 * yp - 2 * ys) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1) / 0.2e1;  //dtheta/dyp
  D(2,1)= -(-pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.1e1 / 0.2e1) - (double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (-2 * zs + 2 * zp) / 0.2e1) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1); //dtheta/dzp
  D(3,1)=(double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (2 * xs - 2 * xp) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1) / 0.2e1;   //dtheta/dxs 
  D(4,1)= (double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (2 * yp - 2 * ys) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1) / 0.2e1;  //dtheta/dys 
  D(5,1)= -(-pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.1e1 / 0.2e1) - (double) (zs - zp) * pow((double) (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp), -0.3e1 / 0.2e1) * (double) (-2 * zs + 2 * zp) / 0.2e1) * pow((double) (1 - (int) pow((double) (zs - zp), (double) 2) / (xs * xs - 2 * xs * xp + xp * xp + ys * ys - 2 * ys * yp + yp * yp + zs * zs - 2 * zs * zp + zp * zp)), -0.1e1 / 0.2e1);   //dtheta/dzs           

  D(6,2)=1;//a1 px
  D(7,3)=1;//a1 py
  D(8,4)=1;//a1 pz
  D(9,5)=1;//a1 m
  D(10,6)=1;//nu px
  D(11,7)=1;//nu py
  D(12,8)=1;//[nu pz

  /*if(debug){
    std::cout << "Derivative for Error propogation " << std::endl; 
    D.Print();
    }*/

  // configure custom weights for determining convergence
  w.ResizeTo(npar);
  for(int i=0;i<w.GetNrows();i++){w(i)=1;}

  if(hasa1 && hasnu){
    // store cov
    cov.ResizeTo(npar,npar);
    incov_temp.Print();
    cov=incov_temp.SimilarityT(D);
    cov_first.ResizeTo(npar,npar);
    cov_first=cov;
    // store par
    par.ResizeTo(npar);
    par=inpar;
    par_first.ResizeTo(npar);
    par_first=par;
    /*if(debug){
      std::cout << "first par " << std::endl; 
      par_first.Print(); 
      std::cout << "first cov " << std::endl;
      cov_first.Print();
      }*/
    field=inStates.front().magneticField();
    //if(debug) std::cout << "TauA1NuNumericalKinematicConstraint::ConfigureIntialState end" << std::endl;
    return true;
  }
  //if(debug) std::cout << "TauA1NuNumericalKinematicConstraint::ConfigureIntialState end" << std::endl;
  return false;
}

std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > TauA1NuNumericalKinematicConstraint::ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint){
  //if(debug){
  // std::cout << "TauA1NuNumericalKinematicConstraint::ConvertStateToParameters start" << std::endl;
    //par.Print();
    //cov.Print();
    //}
  //making refitted states of Kinematic Particles
  std::vector<KinematicState> ns;
  AlgebraicVector7            newParA1;
  AlgebraicSymMatrix          newA1cov(7,0);
  AlgebraicMatrix             ns_cov(inStates.size()*npardim+nposdim,inStates.size()*npardim+nposdim,0);

  // setup a1
  bool hasa1(false);
  bool hasnu(false);
  for(unsigned int s=0; s<inStates.size();s++){
    if(inStates.at(s).particleCharge()!=0){
      hasa1=true;
      AlgebraicVector newPar = asHepVector<7>(inStates.at(s).kinematicParameters().vector());
      for(int i =par_vx;i<=par_m; i++){newParA1(i)=newPar(i+1);}
      // now add updated momentum
      for(int i =par_px; i<=par_m; i++){newParA1(i) = par(i+(a1_px-par_px));}
      // now add cov
      newA1cov = asHepMatrix<7>(inStates.at(s).kinematicParametersError().matrix());
      for(int i = par_vx; i<=par_m; i++){
        for(int j = par_vx; j<=par_m; j++){
	  if(i>=par_px && j>=par_px) newA1cov(i+1,j+1)=cov(i+(a1_px-par_px),j+(a1_px-par_px));
          ns_cov(i+s*npardim+1,j+s*npardim+1)=newA1cov(i+1,j+1);
        }
      }
      KinematicParameters nrPar(newParA1);
      KinematicParametersError nrEr(asSMatrix<7>(newA1cov));
      TrackCharge chl = inStates.at(s).particleCharge();
      KinematicState newState(nrPar,nrEr,chl, field);
      ns.push_back(newState);
    }
  }
  // setup nu
  for(unsigned int s=0; s<inStates.size();s++){
    if(inStates.at(s).particleCharge()==0){
      hasnu=true;
      AlgebraicVector7 newParNu;
      // now add updated momentum
      for(int i =par_vx;i<=par_vz; i++){newParNu(i) = newParA1(i+1);}
      for(int i =par_px; i<=par_pz; i++){newParNu(i) = par(i+(a1_px-par_px)+1);}
      newParNu(par_m)=0;

      // now add cov                                                                                                                                                                                                                        
      AlgebraicSymMatrix newNucov(7,0);

      // copy vertex info
      for(int i = par_vx; i<=par_vz; i++){
	for(int j = par_vx; j<=par_vz; j++){
	  newNucov(i+1,j+1)=newA1cov(i+1,j+1);
	  ns_cov(i+s*npardim+1,j+s*npardim+1)=newNucov(i+1,j+1);
	}
      }
      //copy nu momentum
      for(int i = par_px; i<=par_m; i++){
        for(int j = par_px; j<=par_m; j++){
          if(i>=par_px && j>=par_px && i!=par_m && j!=par_m) newNucov(i+1,j+1)=cov(i+(nu_px-par_px),j+(nu_px-par_px));
	  if(i==par_m && j==par_m) newNucov(i+1,j+1)=1.0E-10;
	    ns_cov(i+s*npardim+1,j+s*npardim+1)=newNucov(i+1,j+1);
        }
      }
      KinematicParameters nrPar(newParNu);
      KinematicParametersError nrEr(asSMatrix<7>(newNucov));
      TrackCharge chl = inStates.at(s).particleCharge();
      KinematicState newState(nrPar,nrEr,chl, field);
      ns.push_back(newState);
    }
  }

  //making resulting vertex
  AlgebraicSymMatrix pCov = newA1cov.sub(1,3);
  float ndf=3.0;
  GlobalPoint vPos (newParA1(par_vx),newParA1(par_vy),newParA1(par_vz));
  VertexState st(vPos,GlobalError( asSMatrix<3>(pCov)));
  KinematicVertexFactory vFactory;
  RefCountedKinematicVertex rVtx = vFactory.vertex(st,chi2,ndf);

  //copy nu momentum                                                                                                                                                                                                                     
  for(int i = par_vx; i<=par_vz; i++){
    for(int j = par_vx; j<=par_vz; j++){
      ns_cov(i+2*npardim+1,j+2*npardim+1)=pCov(i+1,j+1);
    }
  }

  //std::cout << "ns_cov " << ns_cov << std::endl;
  std::pair<std::vector<KinematicState>, AlgebraicMatrix> ns_m(ns,ns_cov);
  //if(debug)std::cout << "TauA1NuNumericalKinematicConstraint::ConvertStateToParameters end" << std::endl;
  return std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix>, RefCountedKinematicVertex >(ns_m,rVtx);
}
