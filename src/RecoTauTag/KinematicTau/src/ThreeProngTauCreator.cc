#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

int ThreeProngTauCreator::create(unsigned int &ambiguity,SelectedKinematicDecay &KFTau){
  //std::cout << "Fit start" << std::endl;
  std::vector<RefCountedKinematicParticle> *pions = new std::vector<RefCountedKinematicParticle>; //3 particles (3pi)
  std::vector<RefCountedKinematicParticle> *neutrinos = new std::vector<RefCountedKinematicParticle>; //1 or 2 particles due to ambiguity (nuGuess1 + nuGuess2)

  std::vector<RefCountedKinematicParticle> *a1 = new std::vector<RefCountedKinematicParticle>;    // a1
  std::vector<RefCountedKinematicParticle> *daughters = new std::vector<RefCountedKinematicParticle>; // nu  + a1

  //std::cout << "Fit A" << std::endl;
  if(!createStartScenario(ambiguity,KFTau, *pions, *neutrinos, *a1)) return 0;
  if (pions->size()!=3 ||(neutrinos->size()!=1 && neutrinos->size()!=2) ){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::create: wrong daughter size. found "<<daughters->size()<<" pis and "<<neutrinos->size()<<" nus. Skip this tauCand";
    return 0;
  }
  //std::cout << "Fit A" << std::endl;
  daughters->push_back(a1->at(0));
  daughters->push_back(neutrinos->at(0));
  //std::cout << "Fit B" << std::endl;
  // daughters->push_back(pions->at(0));
  // pions->push_back(neutrinos->at(0));
  /*
  for(unsigned int i=0;i<daughters->size();i++){
    std::cout << "Current (before) " 
	      << daughters->at(i)->currentState().globalPosition().x() << " " 
	      << daughters->at(i)->currentState().globalPosition().y() << " " 
	      << daughters->at(i)->currentState().globalPosition().z() << " "
	      << daughters->at(i)->currentState().globalMomentum().x() << " "
              << daughters->at(i)->currentState().globalMomentum().y() << " "
              << daughters->at(i)->currentState().globalMomentum().z() << " " <<  daughters->at(i)->currentState().particleCharge() << std::endl; 
  } 
  */
  bool fitWorked=false;
  // fitWorked=kinematicRefit(ambiguity,*pions, modifiedPV_);
  fitWorked=kinematicRefit(ambiguity,*daughters, modifiedPV_);

  /*

  for(unsigned int i=0;i<daughters->size();i++){
    std::cout << "Inital "
              << daughters->at(i)->initialState().globalPosition().x() << " "
              << daughters->at(i)->initialState().globalPosition().y() << " "
              << daughters->at(i)->initialState().globalPosition().z() << " "
              << daughters->at(i)->initialState().globalMomentum().x() << " "
              << daughters->at(i)->initialState().globalMomentum().y() << " "
              << daughters->at(i)->initialState().globalMomentum().z() << " " << daughters->at(i)->currentState().particleCharge() << std::endl;

    std::cout << "Current "
              << daughters->at(i)->currentState().globalPosition().x() << " "
              << daughters->at(i)->currentState().globalPosition().y() << " "
              << daughters->at(i)->currentState().globalPosition().z() << " "
              << daughters->at(i)->currentState().globalMomentum().x() << " "
              << daughters->at(i)->currentState().globalMomentum().y() << " "
              << daughters->at(i)->currentState().globalMomentum().z() << " " << daughters->at(i)->currentState().particleCharge() << std::endl;
  }

  std::cout << "Top Particle "
	    <<  kinTree_->topParticle()->currentState().globalPosition().x() << " "
	    <<  kinTree_->topParticle()->currentState().globalPosition().y() << " "
	    <<  kinTree_->topParticle()->currentState().globalPosition().z() << " "
	    <<  kinTree_->topParticle()->currentState().globalMomentum().x() << " "
	    <<  kinTree_->topParticle()->currentState().globalMomentum().y() << " "
	    <<  kinTree_->topParticle()->currentState().globalMomentum().z() << " " << std::endl;

  TLorentzVector tau_d(kinTree_->topParticle()->currentState().globalMomentum().x(),
		       kinTree_->topParticle()->currentState().globalMomentum().y(),
		       kinTree_->topParticle()->currentState().globalMomentum().z(),
		       sqrt(kinTree_->topParticle()->currentState().globalMomentum().mag2()+pow(kinTree_->topParticle()->currentState().mass(),2.0)));

  if(GenPart.isValid()){                                                                                                                                                                                                                   
    double TauMatchingDR_=0.4;                                                                                                                                                                                                             
    for(reco::GenParticleCollection::const_iterator itr = GenPart->begin(); itr!= GenPart->end(); ++itr){                                                                                                                                  
      if(itr->pdgId()==15){                                                                                                                                                                                                                
	const reco::GenParticle mytau=(*itr);                                                                                                                                                                                              
	TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());                                                                                                                                                     
	if(tau_d.DeltaR(mc)<TauMatchingDR_){                                                                                                                                                                                                
	  TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());                                                                                                                                                   
	  for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){                                                                                                                                                                        
	    const reco::Candidate *dau=(itr)->daughter(i);                                                                                                                                                                                 
	    if(fabs(dau->pdgId())==20213){                                                                                                                                                                                                 
	      std::cout << "Tau Truth " << std::endl;
	      mc.Print();
	      std::cout << "a1 Truth " << std::endl;
	      TLorentzVector a1_d(dau->p4().Px(),dau->p4().Py(),dau->p4().Pz(),dau->p4().E());
	      a1_d.Print();
	    }
	    if(fabs(dau->pdgId())==16){
	      std::cout << "Neutrino Truth " << std::endl;
	      TLorentzVector nu_d(dau->p4().Px(),dau->p4().Py(),dau->p4().Pz(),dau->p4().E());
              nu_d.Print();
	    }
	  }
	}
      }
    }
  } 
  */

  //std::cout << "Fit done" << std::endl;
  delete daughters;
  delete neutrinos;
  delete a1;

  if(fitWorked) return 1;
  else return 0;
}

bool ThreeProngTauCreator::createStartScenario(unsigned int &ambiguity,SelectedKinematicDecay &KFTau, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos, std::vector<RefCountedKinematicParticle> &a1){
  KinematicParticleFactoryFromTransientTrack kinFactory;
  float piMassSigma = sqrt(pow(10.,-12.));//not to small to avoid singularities
  float piChi = 0., piNdf = 0.;//only initial values
  
  // Obtain stored values from Kinematic Tau Candidate
  SecondaryVertexHelper SVH(transientTrackBuilder_,KFTau);
  selectedTracks_=KFTau.InitialTrackTriplet(); //uses orignal triplet tracks (why?)
  std::vector<reco::TransientTrack> transTrkVect=SVH.InitialRefittedTracks();
  TransientVertex secVtx=SVH.InitialSecondaryVertex();
  TLorentzVector lorentzA1=KFTau.Initial_a1_p4();
  double thetaMax=KFTau.InitialThetaMax();
  double theta0=KFTau.InitialThetaGJ();
  TVector3 tauFlghtDir;
  if(ambiguity==SelectedKinematicDecay::PlusSolution || ambiguity==SelectedKinematicDecay::MinusSolution){
    //std::cout << "ambiguity==SelectedKinematicDecay::PlusSolution || ambiguity==SelectedKinematicDecay::MinusSolution " << std::endl;
    tauFlghtDir=KFTau.InitialTauFlightDirectionReFitPvtxToSvtx();
    modifiedPV_=KFTau.InitialPrimaryVertexReFit();
    //tauFlghtDir.Print();
  }
  else{
    //std::cout << "Ambiguity solution" << std::endl;
    tauFlghtDir=KFTau.InitialTauFlightDirectionReFitandRotatedPvtxToSvtx();
    modifiedPV_=KFTau.InitialPrimaryVertexReFitAndRotated();
    //tauFlghtDir.Print();
  }

  //std::cout << "KFTau.InitialTauFlightDirectionReFitPvtxToSvtx then KFTau.InitialTauFlightDirectionReFitandRotatedPvtxToSvtx" << std::endl;
  //KFTau.InitialTauFlightDirectionReFitPvtxToSvtx().Print();
  //KFTau.InitialTauFlightDirectionReFitandRotatedPvtxToSvtx().Print();
  TVector3 startingtauFlghtDir=tauFlghtDir;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now setup Fit parameters
  for(unsigned int i = 0; i!=transTrkVect.size();i++){
    pions.push_back(kinFactory.particle(transTrkVect[i],PMH.Get_piMass(),piChi,piNdf,secVtx.position(),piMassSigma));
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Quality cuts on InitialThetaGJ
  if(theta0>TMath::Pi()/2){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Unrealistic GJ angle ("<<theta0<<"). tauCand skipped!";
    return false;
  }
  
  if ((thetaMax - theta0) < -pow(10.0,-10.0)*thetaMax && thetaMax>0) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: thetaGJ = "<<theta0<<" but replaced by thetaMax = "<<thetaMax;
    theta0 = thetaMax;
  }

  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: rotated PV ("<<modifiedPV_.x()<<","<<modifiedPV_.y()<<","<<modifiedPV_.z()<<"), SV ("<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<"), phi(vtxLink) "<<atan((secVtx.position().y()-modifiedPV_.position().y())/(secVtx.position().x()-modifiedPV_.position().x()))<<", theta "<<theta0;
  
  //std::cout <<"ThreeProngTauCreator::createStartScenario: rotated PV ("<<modifiedPV_.x()<<","<<modifiedPV_.y()<<","<<modifiedPV_.z()<<"), SV ("<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<"), phi(vtxLink) "<<atan((secVtx.position().y()-modifiedPV_.position().y())/(secVtx.position().x()-modifiedPV_.position().x()))<<", theta "<<theta0 << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now create a1 out of pions

  a1 = a1maker(pions,PostFitPions_);
  // std::cout << "a1 from Kalman Fit " << std::endl;
  //lorentzA1.Print();
  if(a1.size()>0)  lorentzA1.SetPxPyPzE(a1.at(0)->currentState().globalMomentum().x(),
					a1.at(0)->currentState().globalMomentum().y(),
					a1.at(0)->currentState().globalMomentum().z(),
					sqrt(a1.at(0)->currentState().globalMomentum().mag2()+pow(a1.at(0)->currentState().mass(),2.0)));
  //std::cout << "a1 from Vertex helper " << std::endl;
  //lorentzA1.Print();

//   std::cout<<"a1 parameters x "<< a1.at(0)->currentState().globalMomentum().x() <<std::endl;
//   std::cout<<"a1 parameters y "<< a1.at(0)->currentState().globalMomentum().y() <<std::endl;
//    std::cout<<"a1 parameters z "<< a1.at(0)->currentState().globalMomentum().z() <<std::endl;
//    std::cout<<"a1 parameters charge "<< a1.at(0)->currentState().particleCharge() <<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now setup the tau and the neutrino
  // clasical configuration
  double tauSolutions = getTauMomentumMagnitudes(ambiguity,lorentzA1.M(),lorentzA1.P(),PMH.Get_tauMass(),thetaMax);
  if (tauSolutions< 0. ) {
    return false;
  } 

  TLorentzVector TauGuessLV;
  TLorentzVector NuGuessLV;
  tauFlghtDir = tauFlghtDir.Unit();
  TauGuessLV.SetXYZM(tauFlghtDir.X()*tauSolutions, tauFlghtDir.Y()*tauSolutions, tauFlghtDir.Z()*tauSolutions, PMH.Get_tauMass());
  //  std::cout << "Default: Px " << TauGuessLV.Px() << " Py " << TauGuessLV.Py() << " Pz " << TauGuessLV.Pz() << " E " << TauGuessLV.E() 
  //    << "ambiguity  " << ambiguity << std::endl;

  // Physical configuration
  if(ambiguity==SelectedKinematicDecay::PlusSolution || ambiguity==SelectedKinematicDecay::MinusSolution){
    if(fabs(lorentzA1.Angle(tauFlghtDir))>=fabs(thetaMax)){
      //case when initial conditions are unphysical reset with 0.XX theta_{GF,max} for starting condition only  
      reco::Vertex tmppVtx=modifiedPV_;
      TransientVertex tmpVtx=secVtx;
      double tmptheta0=theta0;
      VertexRotation vtxC(lorentzA1);
      //vtxC.rotatePV(tmppVtx,tmpVtx,tmptheta0,startingtauFlghtDir,0.975);
    }
    TLorentzVector TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2;
    SolvebyRotation(startingtauFlghtDir,lorentzA1,TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2);
    if(ambiguity==SelectedKinematicDecay::PlusSolution)  TauGuessLV=TauGuessLV1;
    if(ambiguity==SelectedKinematicDecay::MinusSolution) TauGuessLV=TauGuessLV2;
      if(TauGuessLV1.E()!=TauGuessLV2.E()){
	/*	std::cout << "E new: Px " <<  TauGuessLV.Px() << " Py " << TauGuessLV.Py() << " Pz " << TauGuessLV.Pz() << " E " << TauGuessLV.E() 
		  << " ambiguity  " << ambiguity << " angle " << TauGuessLV.Angle(startingtauFlghtDir) 
		  << " thetaGF " << fabs(lorentzA1.Angle(tauFlghtDir)) << " thetaGFmax " << fabs(thetaMax) 
		  << " Dir " <<  startingtauFlghtDir.X() << " " << startingtauFlghtDir.Y() << " " << startingtauFlghtDir.Z() << std::endl;*/
    }
  }
  neutrinos.push_back(unknownNu(TauGuessLV, lorentzA1, secVtx,NuGuessLV));
  //  std::cout << "Sec Vtx " <<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()
  //	    << " Primary " << modifiedPV_.position().x() << " " << modifiedPV_.position().y() << " " << modifiedPV_.position().z() << std::endl;
  KFTau.SetInitialGuess(ambiguity,TauGuessLV,NuGuessLV,startingtauFlghtDir); 
  ///////////////////////////////////////////////////////
  if(neutrinos.size() != 1){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad neutrino size = "<<neutrinos.size();
    return false;
  }
  /*  std::cout << ambiguity << " Tau E " << TauGuessLV.E() << " (" <<  TauGuessLV.Px() << "," <<  TauGuessLV.Py() << "," <<  TauGuessLV.Pz() 
	    << ") Nu E " << NuGuessLV.E()  << " (" << NuGuessLV.Px() << "," << NuGuessLV.Py() << "," << NuGuessLV.Pz() << ")" 
	    << std::endl;*/
  return true;
}

bool ThreeProngTauCreator::kinematicRefit(unsigned int &ambiguity,std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex & primaryVertex){
  if(unfitDaughters.size()!=2){
    edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
    return false;
  }
  // Setup Constraint
  MultiTrackNumericalKinematicConstraint *TauA1NU=new TauA1NuNumericalKinematicConstraint(primaryVertex,PMH.Get_tauMass(),GenPart,1.0,true);
  try{
    kinTree_ = kcvFitter_->fit(unfitDaughters,TauA1NU);//,&vtxGuess);
  }
  catch(VertexException){//("KinematicStatePropagator without material::propagation failed!")
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: VertexException. Skip tau candidate.";
    return false;
  }
  //delete constraint
  delete TauA1NU;
  // Test whether the fit is valid. This is mainly due to unconverged fits.
  if (kinTree_->isValid()) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Valid tree.";
    return true;
  } 
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Warning! Tree is not valid. Skip tauCand.";//DEBUG
  //edm::LogVerbatim("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: ERROR! Tree is not valid. Skip tauCand.";//INFO
  return false;
}

double ThreeProngTauCreator::getTauMomentumMagnitudes(unsigned int& ambiguity,double ma1,double pa1,double M,double theta){
  // return -1. in case of a problem
  double tauMomentumMagnitude = -1.;
  double denominator = (2*pow(ma1,2) + 2*pow(pa1,2)*pow(sin(theta),2));
  if (denominator==0.) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Bad tau magnitude due to zero denominator. Return 0.";
    return tauMomentumMagnitude;
  }
  
  double root = (pow(ma1,2) + pow(pa1,2))*(pow(pow(ma1,2) - pow(M,2),2) -4*pow(M,2)*pow(pa1,2)*pow(sin(theta),2));
  if (root < 0.) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Bad tau magnitude due to negative root. Skip the root part and calculate the solution for a maximal GJ angle. Both solutions will be equal.";
    root = 0.;
  } 
  else root = sqrt(root);
  
  double numerator = (pow(ma1,2) + pow(M,2))*pa1*cos(theta);
  if(ambiguity==SelectedKinematicDecay::ZeroAmbiguitySolution) tauMomentumMagnitude = numerator/denominator;
  if(ambiguity==SelectedKinematicDecay::PlusSolution)          tauMomentumMagnitude = (numerator + root) / denominator;
  if(ambiguity==SelectedKinematicDecay::MinusSolution)         tauMomentumMagnitude = (numerator - root) / denominator;
  
  //catch negatives and 'nan' and negatives.
  if(!(tauMomentumMagnitude>0.)) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Still bad tau magnitude1 of "<<tauMomentumMagnitude <<" replaced by "<<0.;        
    tauMomentumMagnitude = -1.;
  }
  return tauMomentumMagnitude;
}

RefCountedKinematicParticle ThreeProngTauCreator::unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx,TLorentzVector &NuGuessLV){
  NuGuessLV = tauGuess-a1;
  LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::unknownNu: nuGuess (vx, vy, vz, px,py,pz,m) "<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<","<<NuGuessLV.Px()<<","<<NuGuessLV.Py()<<","<<NuGuessLV.Pz()<<","<<NuGuessLV.M()<<", phi: "<<NuGuessLV.Phi();
  if(tauGuess.P()==0.)  NuGuessLV.SetXYZM(1,1,1,0);
  return virtualKinematicParticle(secVtx, NuGuessLV);
}

RefCountedKinematicParticle ThreeProngTauCreator::virtualKinematicParticle(const TransientVertex & vtxGuess, const TLorentzVector & nuGuess){
  const KinematicParameters parameters(AlgebraicVector7(vtxGuess.position().x(),vtxGuess.position().y(),vtxGuess.position().z(),nuGuess.Px(),nuGuess.Py(),nuGuess.Pz(),nuGuess.M()));
  AlgebraicSymMatrix77 Cov;
  AlgebraicSymMatrix33 pvCov=vtxGuess.positionError().matrix();
  for(int i=0; i<MultiTrackNumericalKinematicConstraint::npardim; i++){
    for(int j=0; j<=i; j++){
      if(i<MultiTrackNumericalKinematicConstraint::nposdim) Cov(i,j)=pvCov(i,j);
      else Cov(i,j)=0;
    }
    double par=0;
    if(i==MultiTrackNumericalKinematicConstraint::par_px) par=pow(nuGuess.Px(),2);
    if(i==MultiTrackNumericalKinematicConstraint::par_py) par=pow(nuGuess.Py(),2);
    if(i==MultiTrackNumericalKinematicConstraint::par_pz) par=pow(nuGuess.Pz(),2);
    if(par<25.0) par=25.0;
    Cov(i,i)+=par;
  }

  const KinematicParametersError parametersError(Cov);
  const TrackCharge charge = 0;
  KinematicState kineState(parameters,parametersError,charge,transientTrackBuilder_->field());
  float chiSquared=0.0, degreesOfFr=0.0;
  VirtualKinematicParticleFactory factory;
  return factory.particle(kineState,chiSquared,degreesOfFr,0,0);
}

int ThreeProngTauCreator::ndf() const {
  int freeParameters = 28-15-4-3;
  int constraints = (int)kinTree_->topParticle()->degreesOfFreedom();
  int ndf = constraints - freeParameters;
  //if(ndf != 2) printf("ThreeProngTauCreator::ndf: Warning! Unexpected ndf of %i. Expected 2.\n", ndf );
  return ndf;
}


void ThreeProngTauCreator::quadratic(double &x_plus,double &x_minus,double a, double b, double c){
  double R=b*b-4*a*c;
  if(R<0){R=0;}
  x_minus=(-b-sqrt(R))/(2.0*a);
  x_plus=(-b+sqrt(R))/(2.0*a);
}

void ThreeProngTauCreator::AnalyticESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1){
  double mtau=PMH.Get_tauMass();
  double a=(A1.Pz()*A1.Pz())/(A1.E()*A1.E())-1.0;
  double K=(mtau*mtau-A1.M2()-2.0*A1.Pt()*A1.Pt())/(2.0*A1.E());
  double b=2.0*K*A1.Pz()/A1.E();
  double c=K*K-A1.Pt()*A1.Pt();
  double z1(0),z2(0);
  quadratic(z1,z2,a,b,c);
  nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),z1,sqrt(z1*z1+A1.Pt()*A1.Pt()));
  nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),z2,sqrt(z2*z2+A1.Pt()*A1.Pt()));
}

void ThreeProngTauCreator::NumericalESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1){
  double rmin(-100), rmax(100), step(0.01), mtau2(pow(PMH.Get_tauMass(),2.0)), z1(-999), z2(-999), zmin(-999), min(9999), prev(9999);
  double z=rmin;
  TLorentzVector nu,tau;
  for(int i=0;i<=(int)(rmax-rmin)/step;i++){
    nu.SetPxPyPzE(-A1.Px(),-A1.Py(),z,sqrt(z*z+A1.Pt()*A1.Pt()));
    tau=A1+nu;
    double m2=tau.M2();
    if(m2-mtau2<0 && prev-mtau2>=0) z1=z;
    if(m2-mtau2>0 && prev-mtau2<=0) z2=z;
    if(min>m2){ zmin=z; min=m2;}
    prev=m2;
    z+=step;
  }
  if(z1!=-999 && z2!=-999){
    nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),z1,sqrt(z1*z1+A1.Pt()*A1.Pt()));
    nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),z2,sqrt(z2*z2+A1.Pt()*A1.Pt()));
  }
  else{
    nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),zmin,sqrt(zmin*zmin+A1.Pt()*A1.Pt()));
    nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),zmin,sqrt(zmin*zmin+A1.Pt()*A1.Pt()));
  }
}

void ThreeProngTauCreator::SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau1,TLorentzVector &Tau2,
					   TLorentzVector &nu1,TLorentzVector &nu2){
  TLorentzVector A1rot=A1;
  double phi(TauDir.Phi()),theta(TauDir.Theta());
  A1rot.RotateZ(-phi);
  A1rot.RotateY(-theta);
  /////////////////////////////////////////////////////
  //  NumericalESolver(nu1,nu2,A1rot); // for debugging AnalyticESolver (slow)
  /* std::cout << "NumericalESolver nu" << std::endl;
  nu1.Print();
  nu2.Print();
  AnalyticESolver(nu1,nu2,A1rot);
  std::cout << "AnalyticESolver nu" << std::endl;
  nu1.Print();
  nu2.Print();*/
  AnalyticESolver(nu1,nu2,A1rot);
  /////////////////////////////////////////////////////
  nu1.RotateY(theta);
  nu1.RotateZ(phi);
  Tau1=A1+nu1;
  //
  nu2.RotateY(theta);
  nu2.RotateZ(phi);
  Tau2=A1+nu2;
  /*
  std::cout << "Solution A1" << std::endl;
  A1.Print();
  std::cout << "Solution nu" << std::endl;
  nu1.Print();
  nu2.Print();
  std::cout << "Solution tau" << std::endl;
  Tau1.Print();
  Tau2.Print();
  TauDir.Print();
  std::cout << " " << theta << " " << phi << std::endl;
  */
}

std::vector<RefCountedKinematicParticle> ThreeProngTauCreator::a1maker(std::vector<RefCountedKinematicParticle> &pions,std::vector<RefCountedKinematicParticle> &PostFitPions){
  std::vector<RefCountedKinematicParticle> a1;
  VirtualKinematicParticleFactory KFfactory;
  KinematicParticleVertexFitter kpvFitter;

  RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
  jpTree->movePointerToTheTop();
  // KinematicState kineState = jpTree->currentParticle()->currentState();
  
  const TrackCharge a1Charge = jpTree->currentParticle()->currentState().particleCharge();
  const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
  const KinematicParametersError parametersError = jpTree->currentParticle()->currentState().kinematicParametersError();
  float chi2 = jpTree->currentParticle()->chiSquared();
  float ndf  = jpTree->currentParticle()->degreesOfFreedom();

  KinematicState KinematicStateOfA1(parameters, parametersError, a1Charge, transientTrackBuilder_->field());
  a1.push_back(KFfactory.particle(KinematicStateOfA1, chi2, ndf, 0,0));

  std::vector<RefCountedKinematicParticle> a1daughters = jpTree->daughterParticles();
  for (std::vector<RefCountedKinematicParticle>::iterator iter=a1daughters.begin(); iter!=a1daughters.end(); ++iter) {
    PostFitPions.push_back((*iter));
  }
  return  a1;
}
