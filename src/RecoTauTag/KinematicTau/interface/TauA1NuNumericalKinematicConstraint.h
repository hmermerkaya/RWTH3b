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



class TauA1NuNumericalKinematicConstraint : public MultiTrackNumericalKinematicConstraint{
 public:
  TauA1NuNumericalKinematicConstraint(TVector3 pv,double mtau,edm::Handle<reco::GenParticleCollection> &GenPart_,double weight=1.0,bool debug_=false):
    MultiTrackNumericalKinematicConstraint(weight),
    pv_inital(pv),
    mtau_c(mtau),
    GenPart(GenPart_),
    debug(debug_)
      {
	if(!GenPart.isValid())debug=false;
      }
  virtual ~TauA1NuNumericalKinematicConstraint(){}

  virtual TauA1NuNumericalKinematicConstraint * clone() const {return new TauA1NuNumericalKinematicConstraint(*this);}
  virtual bool   ConfigureIntialState(const std::vector<KinematicState> inStates,const GlobalPoint& inPoint);
  virtual std::pair<std::pair<std::vector<KinematicState>, AlgebraicMatrix >, RefCountedKinematicVertex > ConvertStateToParameters(const std::vector<KinematicState> &inStates,const GlobalPoint& inPoint);
  virtual int numberOfEquations(){return 3;}

 protected:
  virtual AlgebraicVector Value(AlgebraicVector &v);
    
 private:
  TVector3 pv_inital;
  double mtau_c;
  edm::Handle<reco::GenParticleCollection> &GenPart;
  bool debug;
};
#endif
