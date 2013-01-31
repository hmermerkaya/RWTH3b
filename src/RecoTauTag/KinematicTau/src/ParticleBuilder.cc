#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include <TVector3.h>

TrackParticle ParticleBuilder::CreateTrackParticle(const reco::TrackRef &track, edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const GlobalPoint p){
  // Configured for CMSSW Tracks only
  TMatrixT<double> par(TrackParticle::NHelixPar,1);
  TMatrixTSym<double> cov(TrackParticle::NHelixPar);
  double c=1.0;

  for(int i=0;i<TrackParticle::NHelixPar;i++){
    int icmssw=GetCMSSWTrackParIndex(i);
    if(icmssw<0)continue;
    if(i==TrackParticle::kappa){
      if(track->parameter(icmssw)>0){c=1.0;}
      else{c=-1.0;}
      par(i,0)=fabs(track->parameter(icmssw));
    }
    else{par(i,0)=track->parameter(icmssw);}
    for(int j=0;j<TrackParticle::NHelixPar;j++){
      int jcmssw=GetCMSSWTrackParIndex(j);
      if(jcmssw<0) continue;
      cov(i,j)=track->covariance(icmssw,jcmssw);
    }
  }
  ParticleMassHelper PMH;
  return  TrackParticle(par,cov,c,PdtPdgMini::pi0,PMH.Get_piMass(),transTrackBuilder->field()->inInverseGeV(p).z());
}

int ParticleBuilder::GetCMSSWTrackParIndex(int i){
  if(i==TrackParticle::kappa)    return reco::TrackBase::i_qoverp;
  if(i==TrackParticle::lambda)   return reco::TrackBase::i_lambda;
  if(i==TrackParticle::phi)      return reco::TrackBase::i_phi;
  if(i==TrackParticle::dz)       return reco::TrackBase::i_dxy;
  if(i==TrackParticle::dxy)      return reco::TrackBase::i_dsz;
  return -1;
}



reco::Vertex ParticleBuilder::GetVertex(LorentzVectorParticle p){
  TVector3 v=p.Vertex();
  TMatrixTSym<double> vcov=p.VertexCov();
  reco::Vertex::Point vp(v.X(),v.Y(),v.Z());
  reco::Vertex::Error ve;
  for(int i=0;i<vcov.GetNrows();i++){
    for(int j=0;j<=vcov.GetNrows();j++){ve(i,j)=vcov(i,j);}
  }
  return reco::Vertex(vp,ve);
}
