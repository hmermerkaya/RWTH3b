#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryError.h"
#include <TVector3.h>

TrackParticle ParticleBuilder::CreateTrackParticle(const reco::TrackRef &track, edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const GlobalPoint p, bool fromPerigee,bool useTrackHelixPropogation){
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
        cov(i,j)=transTrk.trajectoryStateClosestToPoint(origin).perigeeError().covarianceMatrix()(i,j);
      }
    }
    par(TrackParticle::NHelixPar,0)=transTrackBuilder->field()->inInverseGeV(p).z();
    SFpar=ConvertCMSSWTrackPerigeeToSFTrackPar(par);
    SFcov=ErrorMatrixPropagator::PropogateError(&ParticleBuilder::ConvertCMSSWTrackPerigeeToSFTrackPar,par,cov);
    if(useTrackHelixPropogation){
      /////////////////////////////////////////////////////////////////
      // correct dxy dz neglecting material and radiative corrections
      //std::cout << "Offical CMS dxy " << par(TrackParticle::dxy,0) << " dz " << par(TrackParticle::dz,0) << " " <<  std::endl;
      double x,y,z,dxy,dz,s,kappa,lambda,phi;
      TMatrixT<double>    freehelix(TrackHelixVertexFitter::NFreeTrackPar,1);
      freehelix(TrackHelixVertexFitter::x0,0)=transTrk.initialFreeState().position().x();
      freehelix(TrackHelixVertexFitter::y0,0)=transTrk.initialFreeState().position().y();
      freehelix(TrackHelixVertexFitter::z0,0)=transTrk.initialFreeState().position().z();
      freehelix(TrackHelixVertexFitter::kappa0,0)=SFpar(TrackParticle::kappa,0);
      freehelix(TrackHelixVertexFitter::lambda0,0)=SFpar(TrackParticle::lambda,0);
      freehelix(TrackHelixVertexFitter::phi0,0)=SFpar(TrackParticle::phi,0);
      TrackHelixVertexFitter::Computedxydz(freehelix,0,kappa,lambda,phi,x,y,z,s,dxy,dz);
      SFpar(TrackParticle::dxy,0) = dxy;
      SFpar(TrackParticle::dz,0)  = dz;
      ////////////////////////////////////////////////////////////////
    }
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
  par(TrackParticle::dxy,0)    = -inpar(aTip,0);
  par(TrackParticle::dz,0)     = inpar(aLip,0);
  return par;
}
