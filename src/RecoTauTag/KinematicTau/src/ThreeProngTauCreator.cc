#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

int ThreeProngTauCreator::create(unsigned int &ambiguity,SelectedKinematicDecay &KFTau){
  std::vector<TrackParticle> pions;
  std::vector<LorentzVectorParticle> neutrino;
  std::vector<LorentzVectorParticle> daughters;
  // Helix fit for a1
  if(!FitA1(KFTau)){return 0;}
  A1=mother();
  // Tau fit
  daughters.push_back(A1);
  if(!FitTau(daughters,PV_,ambiguity)){return 1;}
  Tau=mother();
  // Fit sequence complete
  return 2;
}

bool ThreeProngTauCreator::FitTau(std::vector<LorentzVectorParticle>  &unfitDaughters,const reco::Vertex &primaryVertex,unsigned int &ambiguity){
  if(unfitDaughters.size()!=1){
    edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
    return false;
  }
  // Setup Constraint
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
  TauA1NuConstrainedFitter TauA1NU(ambiguity,unfitDaughters.at(0),pv,pvcov);
  TauA1NU.SetMaxDelta(0.01);
  TauA1NU.SetNIterMax(1000);
  TauA1NU.Fit();
  if (TauA1NU.isConverged()) {
    StoreResults(TauA1NU.ChiSquare(),TauA1NU.NDF(),TauA1NU.CSum(),TauA1NU.NIter(),TauA1NU.NConstraints(),TauA1NU.GetReFitDaughters(),TauA1NU.GetMother());
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Valid tree.";
    return true;
  } 
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Warning! Tree is not converged. Skip tauCand.";//DEBUG
  return false;
}

bool ThreeProngTauCreator::FitA1(SelectedKinematicDecay &KFTau){
  if(useTrackHelixFit_){
    std::vector<TrackParticle> pions;
    PV_=KFTau.InitialPrimaryVertexReFit();
    GlobalPoint pvtx(PV_.position().x(),PV_.position().y(),PV_.position().z());
    std::vector<reco::TrackRef> selectedTracks=KFTau.InitialTrackTriplet();
    for(unsigned int i = 0; i!=selectedTracks.size();i++){
      reco::TransientTrack transTrk=transientTrackBuilder_->build(selectedTracks.at(i));
      pions.push_back(ParticleBuilder::CreateTrackParticle(transTrk,transientTrackBuilder_,pvtx,true,true));
    }

    TVector3 pv(PV_.position().x(),PV_.position().y(),PV_.position().z());
    Chi2VertexFitter chi2v(pions,pv);
    chi2v.Fit();
    double c(0); for(unsigned int i=0;i<pions.size();i++){c+=pions.at(i).Charge();}
    int pdgid=fabs(PdtPdgMini::a_1_plus)*c;

    StoreResults(chi2v.ChiSquare(),chi2v.NDF(),0,0,0,chi2v.GetReFitLorentzVectorParticles(),chi2v.GetMother(pdgid));  
  }
  else{
    // Kalman fit
    SecondaryVertexHelper SVH(transientTrackBuilder_,KFTau);
    TransientVertex secVtx=SVH.SecondaryVertex();
    GlobalPoint sv(secVtx.position().x(),secVtx.position().y(),secVtx.position().z());
    std::vector<reco::TransientTrack> transTrkVect=secVtx.refittedTracks();
    // Use CMSSW Kinematic tools to extract a1
    KinematicParticleFactoryFromTransientTrack kinFactory;
    float piMassSigma = sqrt(pow(10.,-12.));//not to small to avoid singularities
    float piChi = 0., piNdf = 0.;//only initial values
    std::vector<RefCountedKinematicParticle> pions;
    for(unsigned int i = 0; i!=transTrkVect.size();i++){
      pions.push_back(kinFactory.particle(transTrkVect[i],PMH.Get_piMass(),piChi,piNdf,secVtx.position(),piMassSigma));
    }
    KinematicParticleVertexFitter kpvFitter;
    RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
    jpTree->movePointerToTheTop();
    const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
    //const KinematicParametersError parametersError = jpTree->currentParticle()->currentState().kinematicParametersError();
    AlgebraicSymMatrix77 cov=jpTree->currentParticle()->currentState().kinematicParametersError().matrix();    
    // get pions
    double c(0); 
    std::vector<reco::Track> Tracks;
    std::vector<LorentzVectorParticle> ReFitPions;
    for(unsigned int i=0;i<transTrkVect.size();i++){
      c+=transTrkVect.at(i).charge();
      reco::Vertex V=secVtx;
      ReFitPions.push_back(ParticleBuilder::CreateLorentzVectorParticle(transTrkVect.at(i),transientTrackBuilder_,V,true,true));      
    }
    // now covert a1 into LorentzVectorParticle
    TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
    TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
    for(int i = 0; i<7; i++){
      a1_par(i,0)=parameters(i);
      for(int j = 0; j<7; j++){a1_cov(i,j)=cov(i,j);}
    }
    LorentzVectorParticle a1(a1_par,a1_cov,fabs(PdtPdgMini::a_1_plus)*c,c,transientTrackBuilder_->field()->inInverseGeV(sv).z());
    // store the results
    StoreResults(secVtx.totalChiSquared(),secVtx.degreesOfFreedom(),0,0,0,ReFitPions,a1);
  }
  return  true;
}

