#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"

int ThreeProngTauCreator::create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks){
  std::vector<RefCountedKinematicParticle> *daughters = new std::vector<RefCountedKinematicParticle>; //3 particles (3pi)
  std::vector<RefCountedKinematicParticle> *neutrinos = new std::vector<RefCountedKinematicParticle>; //1 or 2 particles due to ambiguity (nuGuess1 + nuGuess2)
  std::vector<reco::TrackRef> input = inputTracks;
  
  if(!createStartScenario(input, *daughters, *neutrinos, primaryVertex)) return 0;
  if (daughters->size()!=3 ||(neutrinos->size()!=1 && neutrinos->size()!=2)){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::create: wrong daughter size. found "<<daughters->size()<<" pis and "<<neutrinos->size()<<" nus. Skip this tauCand";
    return 0;
  }
  
  //in this version createStartScenario always rotates up to thetaMax so that there is always only one solution
  daughters->push_back(neutrinos->at(0));
  
  bool fitWorked = kinematicRefit(*daughters, modifiedPV_);//classic: point to rotated PV
  //bool fitWorked = kinematicRefit(*daughters, primaryVertex);//alternative: point to reduced/unrotated PV
  
  delete daughters;
  delete neutrinos;
  
  if(fitWorked) return 1;
  else return 0;
}

bool ThreeProngTauCreator::createStartScenario(std::vector<reco::TrackRef> &input, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos, const reco::Vertex & primaryVertex){
  KinematicParticleFactoryFromTransientTrack kinFactory;
  float piMassSigma = sqrt(pow(10.,-12.));//not to small to avoid singularities
  float piChi = 0., piNdf = 0.;//only initial values
  
  if (input.size()<3) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad track size = "<<input.size();
    return false;
  }
  if (input.size()>3){
    reco::Vertex pvtx=primaryVertex;
    if(!choose3bestTracks(input, pvtx)) return false;
  }
  else{
    if(!sumCharge(input)){
      LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand due to bad charge sum.";
      return false;
    }
    double massA1 = getInvariantMass(input);
    if(massA1 > 2.0 || massA1 < 3*Get_piMass()){//soft upper value
      LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand due to bad a1 mass.";
      return false;
    }
  }
  
  selectedTracks_ = input;
  std::vector<reco::TransientTrack> transTrkVect = convToTransTrck(selectedTracks_);
  if(transTrkVect.size()!=3){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad transient track size = "<<transTrkVect.size();
    return false;
  }
  TransientVertex secVtx;
  if(!checkSecVtx(transTrkVect, secVtx)){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand.";
    return false;
  }
  double massA1 = getInvariantMass(selectedTracks_);
  TLorentzVector lorentzA1 = getSumTLorentzVec(selectedTracks_, massA1);
  VertexRotation vtxC(lorentzA1, 0);
  double thetaMax = vtxC.calcThetaMax();
  if(lorentzA1.M() > 2.0 || lorentzA1.M() < 3*Get_piMass()){//soft upper value due to neutrino resolution
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad a1 mass = "<<lorentzA1.M()<<". Skip tauCand.";
    return false;
  }
  
  double theta0;
  TVector3 tauFlghtDir;
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: initial PV ("<<primaryVertex.x()<<","<<primaryVertex.y()<<","<<primaryVertex.z()<<"), SV ("<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<"), phi(vtxLink) "<<atan((secVtx.position().y()-primaryVertex.position().y())/(secVtx.position().x()-primaryVertex.position().x()));
  modifiedPV_ = primaryVertex;// do not modify the initial vertex, the rotation is only used for a consistent start scenario
  double significance = vtxC.rotatePV(modifiedPV_, secVtx, theta0, tauFlghtDir);//rotation is forced to reach thetaMax
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: rotation significance is "<<significance;
  
  for(unsigned int i = 0; i!=transTrkVect.size();i++){
    pions.push_back(kinFactory.particle(transTrkVect[i],Get_piMass(),piChi,piNdf,secVtx.position(),piMassSigma));
  }
  
  if(theta0>TMath::Pi()/2){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Unrealistic GJ angle ("<<theta0<<"). tauCand skipped!";
    return false;
  }
  
  if ((thetaMax - theta0) < -pow(10.0,-10.0)*thetaMax && thetaMax>0) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: thetaGJ = "<<theta0<<" but replaced by thetaMax = "<<thetaMax;
    theta0 = thetaMax;
  }
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: rotated PV ("<<modifiedPV_.x()<<","<<modifiedPV_.y()<<","<<modifiedPV_.z()<<"), SV ("<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<"), phi(vtxLink) "<<atan((secVtx.position().y()-modifiedPV_.position().y())/(secVtx.position().x()-modifiedPV_.position().x()))<<", theta "<<theta0;
  
  std::pair<double,double> tauSolutions = getTauMomentumMagnitudes(lorentzA1.M(),lorentzA1.P(),Get_tauMass(),theta0);//use a pair to prepare for ambiguities
  if (tauSolutions.first < 0. || tauSolutions.second < 0.) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip invalid tau solutions = "<<tauSolutions.first<<", "<<tauSolutions.second;
    return false;
  } 
  else {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: tau solutions = "<<tauSolutions.first<<", "<<tauSolutions.second;
  }
  
  bool ambiguity = false;
  if(fabs(tauSolutions.first-tauSolutions.second) > pow(10.0,-6)) ambiguity = true;
  
  TLorentzVector tauGuess1;
  tauFlghtDir = tauFlghtDir.Unit();
  tauGuess1.SetXYZM(tauFlghtDir.X()*tauSolutions.first , tauFlghtDir.Y()*tauSolutions.first , tauFlghtDir.Z()*tauSolutions.first , Get_tauMass());
  neutrinos.push_back(unknownNu(tauGuess1, lorentzA1, secVtx));
  if(ambiguity==true){
    TLorentzVector tauGuess2;
    tauGuess2.SetXYZM(tauFlghtDir.X()*tauSolutions.second, tauFlghtDir.Y()*tauSolutions.second, tauFlghtDir.Z()*tauSolutions.second, Get_tauMass());
    neutrinos.push_back(unknownNu(tauGuess2, lorentzA1, secVtx));
  }
  if(neutrinos.size() != 1 && neutrinos.size() != 2){
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad neutrino size = "<<neutrinos.size();
    return false;
  }
  
  return true;
}

bool ThreeProngTauCreator::kinematicRefit(std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex & primaryVertex){
  if(unfitDaughters.size()!=4){
    edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
    return false;
  }
  
  std::vector<MultiTrackKinematicConstraint* > constraintVector;
  MultiTrackKinematicConstraint *tauMass_c = new  MultiTrackMassKinematicConstraint(Get_tauMass(), unfitDaughters.size());
  constraintVector.push_back(tauMass_c);
  GlobalPoint linP(primaryVertex.x(), primaryVertex.y(), primaryVertex.z());
  //	MultiTrackKinematicConstraint *pointing_c = new MultiTrackPointingKinematicConstraint(linP);
  MultiTrackKinematicConstraint *pointing_c = new MultiTrackVertexLinkKinematicConstraint(linP);
  constraintVector.push_back(pointing_c);
  MultiTrackKinematicConstraint *combiC = new CombinedKinematicConstraint(constraintVector);
  
  GlobalPoint vtxGuess = unfitDaughters[3]->currentState().globalPosition();//nu was created at common/corrected vertex of pions
  try{
    kinTree_ = kcvFitter_->fit(unfitDaughters, combiC);
  }
  catch(VertexException){//("KinematicStatePropagator without material::propagation failed!")
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: VertexException. Skip tau candidate.";
    return false;
  }
  
  delete combiC;
  delete pointing_c;
  delete tauMass_c;
  
  // Test whether the fit is valid. This is mainly due to unconverged fits.
  if (kinTree_->isValid()) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Valid tree.";
    return true;
  } 
  else {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Warning! Tree is not valid. Skip tauCand.";//DEBUG
    //edm::LogVerbatim("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: ERROR! Tree is not valid. Skip tauCand.";//INFO
    return false;
  }
}

std::pair<double,double> ThreeProngTauCreator::getTauMomentumMagnitudes(double ma1,double pa1,double M,double theta){
  // return -1. in case of a problem
  double tauMomentumMagnitude1 = -1.;
  double tauMomentumMagnitude2 = -1.;
  
  double denominator = (2*pow(ma1,2) + 2*pow(pa1,2)*pow(sin(theta),2));
  if (denominator==0.) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Bad tau magnitude due to zero denominator. Return 0.";
    return std::make_pair(tauMomentumMagnitude1,tauMomentumMagnitude2);
  }
  
  double root = (pow(ma1,2) + pow(pa1,2))*(pow(pow(ma1,2) - pow(M,2),2) -4*pow(M,2)*pow(pa1,2)*pow(sin(theta),2));
  if (root < 0.) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Bad tau magnitude due to negative root. Skip the root part and calculate the solution for a maximal GJ angle. Both solutions will be equal.";
    root = 0.;
  } 
  else root = sqrt(root);
  
  double numerator = (pow(ma1,2) + pow(M,2))*pa1*cos(theta);
  
  tauMomentumMagnitude1 = (numerator - root) / denominator;
  tauMomentumMagnitude2 = (numerator + root) / denominator;
  
  //catch negatives and 'nan' and negatives.
  if(!(tauMomentumMagnitude1>0.)) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Still bad tau magnitude1 of "<<tauMomentumMagnitude1<<" replaced by "<<0.;        
    tauMomentumMagnitude1 = -1.;
  }
  if(!(tauMomentumMagnitude2>0.)) {
    LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::getTauMomentumMagnitudes: Still bad tau magnitude2 of "<<tauMomentumMagnitude2<<" replaced by "<<0.;        
    tauMomentumMagnitude2 = -1.;
  }
  
  return std::make_pair(tauMomentumMagnitude1,tauMomentumMagnitude2);
}

RefCountedKinematicParticle ThreeProngTauCreator::unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx){
  TLorentzVector nuGuess = tauGuess-a1;
  LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::unknownNu: nuGuess (vx, vy, vz, px,py,pz,m) "<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<","<<nuGuess.Px()<<","<<nuGuess.Py()<<","<<nuGuess.Pz()<<","<<nuGuess.M()<<", phi: "<<nuGuess.Phi();
  if(tauGuess.P()==0.)  nuGuess.SetXYZM(1,1,1,0);
  return virtualKinematicParticle(secVtx, nuGuess);
}

RefCountedKinematicParticle ThreeProngTauCreator::virtualKinematicParticle(const TransientVertex & vtxGuess, const TLorentzVector & nuGuess){
  VirtualKinematicParticleFactory factory;
  //(x,y,z,p_x,p_y,p_z,m)
  const KinematicParameters parameters(AlgebraicVector7(vtxGuess.position().x(),vtxGuess.position().y(),vtxGuess.position().z(),nuGuess.Px(),nuGuess.Py(),nuGuess.Pz(),nuGuess.M()));//Use MET as initial guess for pt?
  ROOT::Math::SVector<double,28> svector28;
  for(unsigned int i=1; i!=22; i++) svector28(i-1) = 0.0;
  for(unsigned int i=22; i!=28; i++) svector28(i-1) = 0.0;//correlation between mass and momentum/vertex
  for(unsigned int n=1; n!=7; n++) svector28(n*(n+1)/2 - 1) = pow(10.,2.);//diagonals, huge error method
  svector28(27) = pow(10.,-12.);//mass error
  //insert 3prong vertex errors
  svector28[ 0] = vtxGuess.positionError().cxx();
  svector28[ 1] = vtxGuess.positionError().cyx();
  svector28[ 2] = vtxGuess.positionError().cyy();
  svector28[ 3] = vtxGuess.positionError().czx();
  svector28[ 4] = vtxGuess.positionError().czy();
  svector28[ 5] = vtxGuess.positionError().czz();
  
  if (std::abs(nuGuess.Px()) >= 5.0) {
    svector28[ 9] = pow(nuGuess.Px(), 2); //assume an error of 100% of the neutrino momentum component
  } 
  else {
    svector28[ 9] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
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
  
  if(ndf != 2) printf("ThreeProngTauCreator::ndf: Warning! Unexpected ndf of %i. Expected 2.\n", ndf );
  
  return ndf;
}

