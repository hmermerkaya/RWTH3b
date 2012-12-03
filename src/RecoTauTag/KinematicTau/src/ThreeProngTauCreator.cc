#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

int ThreeProngTauCreator::create(unsigned int &ambiguity,SelectedKinematicDecay &KFTau){
  std::cout << "Fit start" << std::endl;
  std::vector<RefCountedKinematicParticle> *pions = new std::vector<RefCountedKinematicParticle>; //3 particles (3pi)
  std::vector<RefCountedKinematicParticle> *neutrinos = new std::vector<RefCountedKinematicParticle>; //1 or 2 particles due to ambiguity (nuGuess1 + nuGuess2)

  std::vector<RefCountedKinematicParticle> *a1 = new std::vector<RefCountedKinematicParticle>;    // a1
  std::vector<RefCountedKinematicParticle> *daughters = new std::vector<RefCountedKinematicParticle>; // nu  + a1

  std::cout << "Fit A" << std::endl;
  if(!createStartScenario(ambiguity,KFTau, *pions, *neutrinos, *a1)) return 0;
  if (pions->size()!=3 ||(neutrinos->size()!=1 && neutrinos->size()!=2) ){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::create: wrong daughter size. found "<<daughters->size()<<" pis and "<<neutrinos->size()<<" nus. Skip this tauCand";
    return 0;
  }
  std::cout << "Fit A" << std::endl;
  daughters->push_back(a1->at(0));
  daughters->push_back(neutrinos->at(0));
  std::cout << "Fit B" << std::endl;
  // daughters->push_back(pions->at(0));
  // pions->push_back(neutrinos->at(0));

  for(unsigned int i=0;i<daughters->size();i++){
    std::cout << "Fit C" << std::endl;
        std::cout << "Current " 
	      << daughters->at(i)->currentState().globalPosition().x() << " " 
	      << daughters->at(i)->currentState().globalPosition().y() << " " 
	      << daughters->at(i)->currentState().globalPosition().z() << " "
	      << daughters->at(i)->currentState().globalMomentum().x() << " "
              << daughters->at(i)->currentState().globalMomentum().y() << " "
              << daughters->at(i)->currentState().globalMomentum().z() << " " << std::endl; 
	std::cout << "Fit D" << std::endl;
    }
  bool fitWorked=false;
  // fitWorked=kinematicRefit(ambiguity,*pions, modifiedPV_);
  fitWorked=kinematicRefit(ambiguity,*daughters, modifiedPV_);

  std::cout << "Fit E" << std::endl;
  for(unsigned int i=0;i<daughters->size();i++){
    std::cout << "Fit F" << std::endl;
    std::cout << "Inital "
              << daughters->at(i)->initialState().globalPosition().x() << " "
              << daughters->at(i)->initialState().globalPosition().y() << " "
              << daughters->at(i)->initialState().globalPosition().z() << " "
              << daughters->at(i)->initialState().globalMomentum().x() << " "
              << daughters->at(i)->initialState().globalMomentum().y() << " "
              << daughters->at(i)->initialState().globalMomentum().z() << " " << std::endl;
    std::cout << "Fit G" << std::endl;
    std::cout << "Current "
              << daughters->at(i)->currentState().globalPosition().x() << " "
              << daughters->at(i)->currentState().globalPosition().y() << " "
              << daughters->at(i)->currentState().globalPosition().z() << " "
              << daughters->at(i)->currentState().globalMomentum().x() << " "
              << daughters->at(i)->currentState().globalMomentum().y() << " "
              << daughters->at(i)->currentState().globalMomentum().z() << " " << std::endl;
    std::cout << "Fit H" << std::endl;
    
  }


  std::cout << "Fit done" << std::endl;
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
    tauFlghtDir=KFTau.InitialTauFlightDirectionReFitPvtxToSvtx();
    modifiedPV_=KFTau.InitialPrimaryVertexReFit();
  }
  else{
    tauFlghtDir=KFTau.InitialTauFlightDirectionReFitandRotatedPvtxToSvtx();
    modifiedPV_=KFTau.InitialPrimaryVertexReFitAndRotated();
  }
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
  
  std::cout <<"ThreeProngTauCreator::createStartScenario: rotated PV ("<<modifiedPV_.x()<<","<<modifiedPV_.y()<<","<<modifiedPV_.z()<<"), SV ("<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<"), phi(vtxLink) "<<atan((secVtx.position().y()-modifiedPV_.position().y())/(secVtx.position().x()-modifiedPV_.position().x()))<<", theta "<<theta0 << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now create a1 out of pions

  a1 = a1maker(pions,PostFitPions_);
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
      vtxC.rotatePV(tmppVtx,tmpVtx,tmptheta0,startingtauFlghtDir,0.975);
    }
    TLorentzVector TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2;
    SolvebyRotation(startingtauFlghtDir,lorentzA1,TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2);
    if(ambiguity==SelectedKinematicDecay::PlusSolution)  TauGuessLV=TauGuessLV1;
    if(ambiguity==SelectedKinematicDecay::MinusSolution) TauGuessLV=TauGuessLV2;
      if(TauGuessLV1.E()!=TauGuessLV2.E()){
	std::cout << "E new: Px " <<  TauGuessLV.Px() << " Py " << TauGuessLV.Py() << " Pz " << TauGuessLV.Pz() << " E " << TauGuessLV.E() 
		  << " ambiguity  " << ambiguity << " angle " << TauGuessLV.Angle(startingtauFlghtDir) 
		  << " thetaGF " << fabs(lorentzA1.Angle(tauFlghtDir)) << " thetaGFmax " << fabs(thetaMax) 
		  << " Dir " <<  startingtauFlghtDir.X() << " " << startingtauFlghtDir.Y() << " " << startingtauFlghtDir.Z() << std::endl;
    }
  }
  neutrinos.push_back(unknownNu(TauGuessLV, lorentzA1, secVtx,NuGuessLV));
  //std::cout << "Sec Vtx " <<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()
  //    << " in " << neutrinos.at(0)->currentState().globalPosition().x() << " " << neutrinos.at(0)->currentState().globalPosition().y() << " " << neutrinos.at(0)->currentState().globalPosition().z() << std::endl;
  KFTau.SetInitialGuess(ambiguity,TauGuessLV,NuGuessLV,startingtauFlghtDir); 
  ///////////////////////////////////////////////////////
  if(neutrinos.size() != 1){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad neutrino size = "<<neutrinos.size();
    return false;
  }
  std::cout << ambiguity << " Tau E " << TauGuessLV.E() << " (" <<  TauGuessLV.Px() << "," <<  TauGuessLV.Py() << "," <<  TauGuessLV.Pz() 
	    << ") Nu E " << NuGuessLV.E()  << " (" << NuGuessLV.Px() << "," << NuGuessLV.Py() << "," << NuGuessLV.Pz() << ")" 
	    << std::endl;
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
    std::cout << "C1" << std::endl;
    kinTree_ = kcvFitter_->fit(unfitDaughters,TauA1NU);//,&vtxGuess);
    std::cout << "C2" << std::endl;
  }
  catch(VertexException){//("KinematicStatePropagator without material::propagation failed!")
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: VertexException. Skip tau candidate.";
    return false;
  }
  std::cout << "D" << std::endl;
  //delete constraint
  delete TauA1NU;
  std::cout << "E" << std::endl;
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
  VirtualKinematicParticleFactory factory;
  //(x,y,z,p_x,p_y,p_z,m)
  const KinematicParameters parameters(AlgebraicVector7(vtxGuess.position().x(),vtxGuess.position().y(),vtxGuess.position().z(),nuGuess.Px(),nuGuess.Py(),nuGuess.Pz(),nuGuess.M()));
  ROOT::Math::SVector<double,28> svector28;
  for(unsigned int i=1; i!=22; i++) svector28(i-1) = 0.0;
  for(unsigned int i=22; i!=28; i++) svector28(i-1) = 0.0; //correlation between mass and momentum/vertex
  for(unsigned int n=1; n!=7; n++) svector28(n*(n+1)/2 - 1) = pow(10.,2.);//diagonals, huge error method
  svector28(27) = pow(10.,-12.);//mass error
  //insert 3prong vertex errors
  svector28[0] = vtxGuess.positionError().cxx()*100;
  svector28[1] = vtxGuess.positionError().cyx()*100;
  svector28[2] = vtxGuess.positionError().cyy()*100;
  svector28[3] = vtxGuess.positionError().czx()*100;
  svector28[4] = vtxGuess.positionError().czy()*100;
  svector28[5] = vtxGuess.positionError().czz()*100;
  
  if (std::abs(nuGuess.Px()) >= 5.0) {
    svector28[9] = pow(nuGuess.Px(), 2); //assume an error of 100% of the neutrino momentum component
  } 
  else {
    svector28[9] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
  }
  if (std::abs(nuGuess.Py()) >= 5.0) {
    svector28[14] = pow(nuGuess.Py(), 2); //assume an error of 100% of the neutrino momentum component
  } 
  else {
    svector28[14] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
  }
  if (std::abs(nuGuess.Pz()) >= 5.0) {
    svector28[20] = pow(nuGuess.Pz(), 2); //assume an error of 100% of the neutrino momentum component
  } 
  else {
    svector28[20] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
  }
  
  ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > matrix(svector28);
  const KinematicParametersError parametersError(matrix);
  const TrackCharge charge = 0;
  KinematicState kineState(parameters, parametersError, charge, transientTrackBuilder_->field());
  float chiSquared=0.0, degreesOfFr=0.0;
  return factory.particle(kineState, chiSquared, degreesOfFr, 0,0);
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

void ThreeProngTauCreator::ESolver(double &Pnuz1,double &Pnuz2,double Ea1,double ma1, double Pz, double Pt){
  double mtau=PMH.Get_tauMass();
  double a=(Pz*Pz)/(Ea1*Ea1)-1.0;
  double K=(mtau*mtau-ma1*ma1-2.0*Pt*Pt)/(2.0*Ea1);
  double b=2.0*K*Pz/Ea1;
  double c=K*K-Pt*Pt;
  quadratic(Pnuz1,Pnuz2,a,b,c);
}


void ThreeProngTauCreator::SolvebyRotation(TVector3 TauDir,TLorentzVector a1, TLorentzVector &Tau1,TLorentzVector &Tau2,
					   TLorentzVector &nu1,TLorentzVector &nu2){
  TLorentzVector a1rot=a1;
  double phi(TauDir.Phi()),theta(TauDir.Theta());
  a1rot.RotateZ(-phi);
  a1rot.RotateY(-theta);
  double NuSolution1(0), NuSolution2(0), Ea1(a1.E()),Pz(a1rot.Pz()),Pt(a1rot.Pt()),ma1(a1.M());
  ESolver(NuSolution1,NuSolution2,Ea1,ma1,Pz,Pt);
  TLorentzVector Neutrino1(-a1rot.Px(),-a1rot.Py(),NuSolution1,sqrt(NuSolution1*NuSolution1+Pt*Pt));
  Neutrino1.RotateY(theta);
  Neutrino1.RotateZ(phi);
  Tau1=a1+Neutrino1;
  nu1=Neutrino1;
  TLorentzVector Neutrino2(-a1rot.Px(),-a1rot.Py(),NuSolution2,sqrt(NuSolution2*NuSolution2+Pt*Pt));
  Neutrino2.RotateY(theta);
  Neutrino2.RotateZ(phi);
  Tau2=a1+Neutrino2;
  a1rot.RotateY(theta);
  a1rot.RotateZ(phi);
  nu2=Neutrino2;
}

std::vector<RefCountedKinematicParticle> ThreeProngTauCreator::a1maker(std::vector<RefCountedKinematicParticle> &pions,std::vector<RefCountedKinematicParticle> &PostFitPions ){

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

  std::vector<RefCountedKinematicParticle> daughters = jpTree->daughterParticles();
  for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter) {
    PostFitPions.push_back((*iter));
  }

  return  a1;
}
