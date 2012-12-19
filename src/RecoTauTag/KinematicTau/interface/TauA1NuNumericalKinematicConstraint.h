#ifndef TauA1NuNumericalKinematicConstraint_H
#define TauA1NuNumericalKinematicConstraint_H

#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "TVector3.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "RecoTauTag/KinematicTau/interface/ThreeProngTauSolver.h"


class TauA1NuNumericalKinematicConstraint : public MultiTrackNumericalKinematicConstraint,  public ThreeProngTauSolver {
 public:
  TauA1NuNumericalKinematicConstraint(unsigned int &ambiguity_,const reco::Vertex &primaryVertex,double mtau,edm::Handle<reco::GenParticleCollection> &GenPart_,double weight=1.0,bool debug_=false);
  virtual ~TauA1NuNumericalKinematicConstraint(){}

  enum Pars{tau_phi=0,tau_theta,a1_px,a1_py,a1_pz,a1_m,nu_px,nu_py,nu_pz,npar,norigpar=13};

  virtual TauA1NuNumericalKinematicConstraint * clone() const {return new TauA1NuNumericalKinematicConstraint(*this);}
  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint);
  virtual int numberOfEquations(){return 3;}

 protected:
  virtual TVectorD Value(TVectorD &v);
    
 private:
  reco::Vertex pv_inital;
  TVector3 TauDir,sv,pv;
  double mtau_c;
  edm::Handle<reco::GenParticleCollection> &GenPart;
  unsigned int ambiguity;

};
#endif
