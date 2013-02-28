#ifndef FitSequencer_h
#define FitSequencer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TString.h"

class FitSequencer {
public:
  FitSequencer(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder,edm::Handle<reco::GenParticleCollection> &GenPart);
  FitSequencer(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const edm::ParameterSet& cfg,edm::Handle<reco::GenParticleCollection> &GenPart);
  virtual ~FitSequencer(){};
  
  virtual int create(unsigned int& ambiguity,SelectedKinematicDecay &KFTau) = 0;

  virtual void  GetSelectedKinematicParticleList(unsigned int &ambiguity,SelectedKinematicParticleCollection &refitTauDecay);
  virtual float ndf(){ return NDF_.at(0);}//double x=0; for(unsigned int i=0;i<NDF_.size() ;i++)  x+=NDF_.at(i);   return x;}
  virtual float chi2(){ return Chi2_.at(0);}//double x=0; for(unsigned int i=0;i<Chi2_.size();i++)  x+=Chi2_.at(i);  return x;}
  virtual float CSum(){ double x=0; for(unsigned int i=0;i<cSum_.size();i++)  x+=cSum_.at(i);  return x;}
  virtual float Niter(){double x=0; for(unsigned int i=0;i<Niter_.size();i++) x+=Niter_.at(i); return x;}
  virtual float NConstraints(){double x=0; for(unsigned int i=0;i<NConstraints_.size();i++) x+=NConstraints_.at(i); return x;}
  virtual LorentzVectorParticle mother(int i=-1);
  virtual std::vector<LorentzVectorParticle> daugthers(int i=-1);
  virtual std::vector<LorentzVectorParticle> chargedDaughters(int i=-1);
  virtual std::vector<LorentzVectorParticle> neutralDaughters(int i=-1);
  virtual reco::PFTau getPFTau();

protected:
  void StoreResults(float chi2,float ndf,float csum,float niter,float nconst,std::vector<LorentzVectorParticle> fitdaughters,LorentzVectorParticle fitmother);

  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder_;
  edm::Handle<reco::GenParticleCollection> &GenPart_;
  edm::ParameterSet cfg;
  LorentzVectorParticle A1,Tau;
  reco::Vertex PV_;
  int status;

 private:
  std::vector<float> Chi2_,NDF_,cSum_,Niter_,NConstraints_;
  std::vector<LorentzVectorParticle> FitMother_;
  std::vector<std::vector<LorentzVectorParticle> > FitDaughters_;
  std::vector<std::vector<LorentzVectorParticle> > unFitDaughters_;
  std::vector<std::vector<TString> > DaughtersNames_;
  std::vector<TString> MothersNames_;
};

#endif
