#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"

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
  // GlobalPoint pv(PV_.position().x(),PV_.position().y(),PV_.position().z());
  std::vector<reco::TrackRef> selectedTracks=KFTau.InitialTrackTriplet();
  ////////// debug
  for(unsigned int i = 0; i!=selectedTracks.size();i++){
    std::cout << "selected tracks px " << selectedTracks.at(i)->px() << " py " << selectedTracks.at(i)->py() << " pz " << selectedTracks.at(i)->pz() << " vx "
              << selectedTracks.at(i)->vx() << " vy " <<selectedTracks.at(i)->vy() << " vz " << selectedTracks.at(i)->vz() << std::endl;

    std::cout << "qoverp " << selectedTracks.at(i)->qoverp() << " theta " << selectedTracks.at(i)->theta() << " kappa" <<  selectedTracks.at(i)->qoverp()/fabs(cos(selectedTracks.at(i)->theta())) 
	      << " d0 " << selectedTracks.at(i)->dxy() << " dsz " << selectedTracks.at(i)->dsz() <<  " dz " << selectedTracks.at(i)->dz() << " cos(lambda)" << cos(selectedTracks.at(i)->lambda()) << " " 
	      << selectedTracks.at(i)->pt()/selectedTracks.at(i)->p() <<   std::endl;  
  }
  SecondaryVertexHelper SVH(transientTrackBuilder_,KFTau);
  TransientVertex secVtx=SVH.InitialSecondaryVertex();
  PV_=secVtx;
  GlobalPoint pv(0.0,0.0,0.0);
  GlobalPoint sv(secVtx.position().x(),secVtx.position().y(),secVtx.position().z());
  std::vector<reco::TransientTrack> transTrkVect=SVH.InitialRefittedTracks();
  GlobalPoint cptpv;
  GlobalPoint cptsv;
  for(unsigned int i = 0; i!=transTrkVect.size();i++){
    std::cout << "transTrkVect tracks px " << transTrkVect.at(i).track().px() << " py " << transTrkVect.at(i).track().py() << " pz " << transTrkVect.at(i).track().pz() << " vx"
              << transTrkVect.at(i).track().vx() << " vy" <<transTrkVect.at(i).track().vy() << " vz " << transTrkVect.at(i).track().vz() << std::endl;
    std::cout << "vx " <<transTrkVect.at(i).trajectoryStateClosestToPoint(pv).position().x() << " vy " << transTrkVect.at(i).trajectoryStateClosestToPoint(pv).position().y() << " vz " << transTrkVect.at(i).trajectoryStateClosestToPoint(pv).position().z() << std::endl;
    cptpv=transTrkVect.at(i).trajectoryStateClosestToPoint(pv).position();
    cptsv=transTrkVect.at(i).trajectoryStateClosestToPoint(sv).position();
    std::cout << "linear estimate of s: " << sqrt(pow(cptpv.x()-cptsv.x(),2.0)+pow(cptpv.y()-cptsv.y(),2.0)+pow(cptpv.z()-cptsv.z(),2.0)) << std::endl;
    std::cout << "pv x0 " << cptpv.x() << " pv y0 " << cptpv.y() << " pv z0 " << cptpv.z() << "pv (x^2+y^2)^1/2 " << sqrt(cptpv.x()*cptpv.x()+cptpv.y()*cptpv.y())  << std::endl;
    std::cout << "sv x0 " << cptsv.x() << " sv y0 " << cptsv.y() << " sv z0 " << cptsv.z() <<  std::endl;
  }
  std::cout << "Kalman Fit Vertex x " << secVtx.position().x() << " y " << secVtx.position().y() << " z " << secVtx.position().z() << std::endl; 
  /////////////// end debug
  for(unsigned int i = 0; i!=selectedTracks.size();i++){
    pions.push_back(ParticleBuilder::CreateTrackParticle(selectedTracks.at(i),transientTrackBuilder_,cptpv));
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
  TLorentzVector LV=chi2v.GetMother(pdgid).LV();
  for(reco::GenParticleCollection::const_iterator itr = GenPart_->begin(); itr!= GenPart_->end(); ++itr){
    if(fabs(itr->pdgId())==15){
      const reco::GenParticle mytau=(*itr);
      for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){
	const reco::Candidate *dau=(itr)->daughter(i);
	if(fabs(dau->pdgId())==20213){
	  TLorentzVector a1(dau->p4().Px(),dau->p4().Py(),dau->p4().Pz(),dau->p4().E());
	  //if(a1.DeltaR(LV)<0.4){
	    // Print info from matching tau
	    std::cout << "A1 Truth Px " << dau->p4().Px() << " Py " << dau->p4().Py() << " Pz " << dau->p4().Pz() << " E " << dau->p4().E() << " M " << dau->p4().M() <<  std::endl;
	    std::cout << "Tau vx "      << (itr)->vx() << " vy " << (itr)->vy() << " vz " << (itr)->vz() << std::endl;
	    std::cout << "A1 vx "       << dau->vx()  << " vy " << dau->vy() << " vz " << dau->vz() << std::endl;
	    //}
	}
      }
    }
  }
  std::cout << "Chi2Vertex Fit vx " << chi2v.GetVertex().X() << " vy " << chi2v.GetVertex().Y() << " " << chi2v.GetVertex().Z() << " " << std::endl; 
  std::cout << "Chi2Vertex a1 Px  " << LV.Px() << " py " << LV.Py() << " pz " << LV.Pz() << std::endl; 
  StoreResults(chi2v.ChiSquare(),chi2v.NDF(),0,0,0,chi2v.GetReFitLorentzVectorParticles(),chi2v.GetMother(pdgid));
  return  true;
}

