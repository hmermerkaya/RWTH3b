#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryError.h"
#include <TVector3.h>

TrackParticle ParticleBuilder::CreateTrackParticle(const reco::TrackRef &track, edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const GlobalPoint p, bool fromPerigee){
  // Configured for CMSSW Tracks only
  reco::TransientTrack transTrk=transTrackBuilder->build(track); 
  TMatrixT<double>    par(TrackParticle::NHelixPar+1,1);
  TMatrixTSym<double> cov(TrackParticle::NHelixPar+1);
  TMatrixT<double>    SFpar(TrackParticle::NHelixPar,1);
  TMatrixTSym<double> SFcov(TrackParticle::NHelixPar);
  if(!fromPerigee){
    for(int i=0;i<TrackParticle::NHelixPar;i++){
      par(i,0)=track->parameter(i);
      for(int j=0;j<TrackParticle::NHelixPar;j++){
	cov(i,j)=track->covariance(i,j);
      }
    }
    par(TrackParticle::NHelixPar,0)=transTrackBuilder->field()->inInverseGeV(p).z();
    SFpar=ConvertCMSSWTrackParToSFTrackPar(par);
    SFcov=ErrorMatrixPropagator::PropogateError(&ParticleBuilder::ConvertCMSSWTrackParToSFTrackPar,par,cov);
  }
  else{
    reco::TransientTrack transTrk=transTrackBuilder->build(track);
    GlobalPoint origin(0.0,0.0,0.0);
    for(int i=0;i<TrackParticle::NHelixPar;i++){
      par(i,0)=transTrk.trajectoryStateClosestToPoint(origin).perigeeParameters().vector()(i);
      for(int j=0;j<TrackParticle::NHelixPar;j++){
        cov(i,j)=transTrk.trajectoryStateClosestToPoint(origin).perigeeParameters().vector()(i);
      }
    }
    par(TrackParticle::NHelixPar,0)=transTrackBuilder->field()->inInverseGeV(p).z();
    SFpar=ConvertCMSSWTrackPerigeeToSFTrackPar(par);
    SFcov=ErrorMatrixPropagator::PropogateError(&ParticleBuilder::ConvertCMSSWTrackPerigeeToSFTrackPar,par,cov);
  
  GlobalPoint pv=transTrk.trajectoryStateClosestToPoint(origin).position();
  std::cout << "BField Z pv " << transTrackBuilder->field()->inInverseGeV(pv).z()  << " Bz/B "<< transTrackBuilder->field()->inInverseGeV(pv).z()/transTrackBuilder->field()->inInverseGeV(pv).mag() << std::endl;
  std::cout << "BField Z sv " << transTrackBuilder->field()->inInverseGeV(p).z()  << " Bz/B "<< transTrackBuilder->field()->inInverseGeV(p).z()/transTrackBuilder->field()->inInverseGeV(p).mag() << std::endl;
  }
  ParticleMassHelper PMH;
  double c=track->charge();
  return  TrackParticle(SFpar,SFcov,abs(PdtPdgMini::pi_plus)*c,PMH.Get_piMass(),c,transTrackBuilder->field()->inInverseGeV(p).z());
}

reco::Vertex ParticleBuilder::GetVertex(LorentzVectorParticle p){
  TVector3 v=p.Vertex();
  TMatrixTSym<double> vcov=p.VertexCov();
  reco::Vertex::Point vp(v.X(),v.Y(),v.Z());
  reco::Vertex::Error ve;
  for(int i=0;i<vcov.GetNrows();i++){
    for(int j=0;j<vcov.GetNrows();j++){ve(i,j)=vcov(i,j);}
  }
  return reco::Vertex(vp,ve);
}


TMatrixT<double> ParticleBuilder::ConvertCMSSWTrackParToSFTrackPar(TMatrixT<double> &inpar){
  TMatrixT<double> par(TrackParticle::NHelixPar,1);
  par(TrackParticle::kappa,0)  = -1.0*inpar(TrackParticle::NHelixPar,0)*inpar(reco::TrackBase::i_qoverp,0)/fabs(cos(inpar(reco::TrackBase::i_lambda,0)));
  par(TrackParticle::lambda,0) = inpar(reco::TrackBase::i_lambda,0);
  par(TrackParticle::phi,0)    = inpar(reco::TrackBase::i_phi,0);
  par(TrackParticle::dz,0)     = inpar(reco::TrackBase::i_dsz,0)/fabs(cos(inpar(reco::TrackBase::i_lambda,0)));
  par(TrackParticle::dxy,0)    = inpar(reco::TrackBase::i_dxy,0);
  return par;
}


TMatrixT<double> ParticleBuilder::ConvertCMSSWTrackPerigeeToSFTrackPar(TMatrixT<double> &inpar){
  TMatrixT<double> par(TrackParticle::NHelixPar,1);
  par(TrackParticle::kappa,0)  = inpar(aCurv,0); 
  par(TrackParticle::lambda,0) = TMath::Pi()/2-inpar(aTheta,0); 
  par(TrackParticle::phi,0)    = inpar(aPhi,0);
  par(TrackParticle::dz,0)     = inpar(aTip,0);
  par(TrackParticle::dxy,0)    = inpar(aLip,0);
  return par;
}
