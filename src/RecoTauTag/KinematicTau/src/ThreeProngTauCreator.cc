#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"

int ThreeProngTauCreator::create(unsigned int &ambiguity,SelectedKinematicDecay &KFTau){
  status=0;
  std::cout << "Fit start" << std::endl;
  std::vector<TrackParticle> pions;
  std::vector<LorentzVectorParticle> neutrino;
  std::vector<LorentzVectorParticle> daughters;
  // Helix fit for a1
  ConfigurePions(KFTau,pions);
  std::cout << "configure pions" << std::endl;
  if(pions.size()!=3)return status;
  std::cout << "pions ok" << std::endl;
  if(!FitA1(pions,PV_)){ std::cout << "status failed" <<  std::endl; return status;}
  std::cout << "vertex and a1 Fit" << std::endl;
  A1=mother();
  std::cout << "have a1 mother" << std::endl;
  // Tau fit

  ConfigureNeutrino(KFTau,ambiguity,A1,neutrino);
  daughters.push_back(A1);
  daughters.push_back(neutrino.at(0));
  std::cout << "have tau daughters" << std::endl;
  if(!FitTau(daughters,PV_,ambiguity)) return status;
  std::cout << " have tau Fit" << std::endl;
  Tau=mother();
  std::cout << " have tau mother" << std::endl;
  // Fit sequence complete
  status=1;
  return status;
}

void ThreeProngTauCreator::ConfigurePions(SelectedKinematicDecay &KFTau, std::vector<TrackParticle> &pions){
  PV_=KFTau.InitialPrimaryVertexReFit();
  std::vector<reco::TrackRef> selectedTracks=KFTau.InitialTrackTriplet();
  GlobalPoint pv(PV_.position().x(),PV_.position().y(),PV_.position().z());
  for(unsigned int i = 0; i!=selectedTracks.size();i++){
    pions.push_back(ParticleBuilder::CreateTrackParticle(selectedTracks.at(i),transientTrackBuilder_,pv));
  }
}

void ThreeProngTauCreator::ConfigureNeutrino(SelectedKinematicDecay &KFTau,int ambiguity,LorentzVectorParticle &a1,std::vector<LorentzVectorParticle> &neutrinos){
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now setup the tau and the neutrino
  TLorentzVector lorentzA1=a1.LV();
  TVector3 sv=a1.Vertex();
  TVector3 pv(PV_.position().x(),PV_.position().y(),PV_.position().z());
  TVector3 tauFlghtDir=sv-pv;
  TLorentzVector TauGuessLV;
  TLorentzVector NuGuessLV;

  // Physical configuration
  TVector3 startingtauFlghtDir=tauFlghtDir.Unit();
  TLorentzVector TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2;
  MultiProngTauSolver Tausolver(PMH.Get_tauMass());
  Tausolver.SolvebyRotation(startingtauFlghtDir,lorentzA1,TauGuessLV1,TauGuessLV2,NuGuessLV1,NuGuessLV2);
  if(ambiguity==SelectedKinematicDecay::PlusSolution)  TauGuessLV=TauGuessLV1;
  if(ambiguity==SelectedKinematicDecay::MinusSolution) TauGuessLV=TauGuessLV2;
  NuGuessLV = TauGuessLV-lorentzA1;
  neutrinos.push_back(EstimateNu(a1,NuGuessLV));
  KFTau.SetInitialGuess(ambiguity,TauGuessLV,NuGuessLV,startingtauFlghtDir); 
}

LorentzVectorParticle ThreeProngTauCreator::EstimateNu(LorentzVectorParticle &a1,TLorentzVector &nuGuess){
  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
  par(LorentzVectorParticle::vx,0)=a1.Parameter(LorentzVectorParticle::vx);
  par(LorentzVectorParticle::vy,0)=a1.Parameter(LorentzVectorParticle::vy);
  par(LorentzVectorParticle::vz,0)=a1.Parameter(LorentzVectorParticle::vz);
  par(LorentzVectorParticle::px,0)=nuGuess.Px();
  par(LorentzVectorParticle::py,0)=nuGuess.Py();
  par(LorentzVectorParticle::pz,0)=nuGuess.Pz();
  par(LorentzVectorParticle::m,0) =nuGuess.M();
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);
  TMatrixTSym<double> pvCov=a1.VertexCov();
  for(int i=0; i<LorentzVectorParticle::NLorentzandVertexPar; i++){
    for(int j=0; j<=i; j++){
      if(i<LorentzVectorParticle::NVertex) Cov(i,j)=pvCov(i,j);
      else Cov(i,j)=0;
    }
    double v=0;
    if(i==LorentzVectorParticle::px || i==LorentzVectorParticle::py || i==LorentzVectorParticle::pz) v=par(i,0)*par(i,0);
    if(v<25.0) v=25.0;
    Cov(i,i)+=v;
  }
  return LorentzVectorParticle(par,Cov,PdtPdgMini::nu_tau,0,a1.BField());
}

bool ThreeProngTauCreator::FitTau(std::vector<LorentzVectorParticle>  &unfitDaughters,const reco::Vertex &primaryVertex,
				  unsigned int &ambiguity){
  if(unfitDaughters.size()!=2){
    edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
    return false;
  }
  // Setup Constraint
  std::cout << "PV " << primaryVertex.position().x() << " " << primaryVertex.position().y() << " " << primaryVertex.position().z() << std::endl;
  TVector3 pv(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);
  math::Error<LorentzVectorParticle::NVertex>::type pvCov;
  primaryVertex.fill(pvCov);
  for(int i = 0; i<LorentzVectorParticle::NVertex; i++){
    for(int j = 0; j<LorentzVectorParticle::NVertex; j++){
      pvcov(i,j)=pvCov(i,j);
      pvcov(j,i)=pvCov(i,j);
    }
  }
  TauA1NuConstrainedFitter TauA1NU(ambiguity,unfitDaughters,pv,pvcov,PMH.Get_tauMass());
  TauA1NU.Fit();
  StoreResults(TauA1NU.ChiSquare(),TauA1NU.NDF(),TauA1NU.CSum(),TauA1NU.NIter(),TauA1NU.NConstraints(),TauA1NU.GetReFitDaughters(),TauA1NU.GetMother());
  if (TauA1NU.isConverged()) {
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Valid tree.";
    return true;
  } 
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Warning! Tree is not converged. Skip tauCand.";//DEBUG
  return false;
}

bool ThreeProngTauCreator::FitA1(std::vector<TrackParticle> &pions,const reco::Vertex & primaryVertex){
  TVector3 pv(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  Chi2VertexFitter chi2v(pions,pv);
  chi2v.Fit();
  double c(0); for(unsigned int i=0;i<pions.size();i++){c+=pions.at(i).Charge();}
  int pdgid=fabs(PdtPdgMini::a_1_plus)*c;
  StoreResults(chi2v.ChiSquare(),chi2v.NDF(),0,0,0,chi2v.GetReFitLorentzVectorParticles(),chi2v.GetMother(pdgid));
  return  true;
}

