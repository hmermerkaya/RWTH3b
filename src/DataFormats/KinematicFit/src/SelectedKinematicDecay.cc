#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

SelectedKinematicDecay::SelectedKinematicDecay() {
  SetInitialProperties(Undefined,reco::PFTauRef(),std::vector<reco::TrackRef>(),reco::Vertex(),"",0);
  SetKinematicFitProperties(SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0,std::map<std::string, bool>());
  SetMissingQualityCriteria(-1,-1,-1,-1);
  SetInitalVertexProperties(reco::Vertex(),reco::Vertex(),std::vector<reco::TransientTrack>(),TransientVertex());  
  SetInitalKinematics(TVector3(),std::vector<TLorentzVector>(),TLorentzVector(),0.0,0.0);
}

SelectedKinematicDecay::SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRefOrig, std::vector<reco::TrackRef> &TrackTriplet,
					       const reco::Vertex &primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx){
  SetInitialProperties(tauDecayMode,tauRefOrig,TrackTriplet,primaryVertex,primVtxReFitTag,nTauPerVtx);
  SetKinematicFitProperties(SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0,std::map<std::string, bool>());
  SetMissingQualityCriteria(-1,-1,-1,-1);
  SetInitalVertexProperties(reco::Vertex(),reco::Vertex(),std::vector<reco::TransientTrack>(),TransientVertex());
  SetInitalKinematics(TVector3(),std::vector<TLorentzVector>(),TLorentzVector(),0.0,0.0);
}

SelectedKinematicDecay::SelectedKinematicDecay(const SelectedKinematicParticleCollection & particles,
					       const int iterations, const int maxiterations, const float csum,
					       const float mincsum, const int constraints, const int ndf, const float chi2,
					       const reco::PFTauRef & tauRef, const std::map<std::string, bool> & discriminators){
  SetInitialProperties(Undefined,tauRef,std::vector<reco::TrackRef>(),reco::Vertex()," ",0);
  SetKinematicFitProperties(particles,iterations,maxiterations,csum,mincsum,constraints,ndf,chi2,discriminators);
  SetMissingQualityCriteria(-1,-1,-1,-1);
  SetInitalVertexProperties(reco::Vertex(),reco::Vertex(),std::vector<reco::TransientTrack>(),TransientVertex());
  SetInitalKinematics(TVector3(),std::vector<TLorentzVector>(),TLorentzVector(),0.0,0.0);
}

SelectedKinematicDecay::~SelectedKinematicDecay(){

}

void SelectedKinematicDecay::SetInitialProperties(unsigned int tauDecayMode,reco::PFTauRef tauRefOrig,
						  std::vector<reco::TrackRef> TrackTriplet, const reco::Vertex primaryVertex,
						  std::string primVtxReFitTag, unsigned int nTauPerVtx){
  tauDecayMode_=tauDecayMode;
  nTauPerVtx_=nTauPerVtx;
  PFTauRefOrig_=tauRefOrig;
  primVtxReFitTag_=primVtxReFitTag;

  initalTrackTriplet_=TrackTriplet;
  initalPrimVtx_=primaryVertex;
}

void SelectedKinematicDecay::SetKinematicFitProperties(const SelectedKinematicParticleCollection particles,
						       const int iterations, const int maxiterations, const float csum,
						       const float mincsum, const int constraints, const int ndf,
						       const float chi2, const std::map<std::string, bool> discriminators){
  particles_ = particles;
  iterations_ = iterations;
  maxiterations_ = maxiterations;
  csum_ = csum;
  mincsum_ = mincsum;
  constraints_ = constraints;
  ndf_ = ndf;
  chi2_ = chi2;
  discriminators_ = discriminators;
}

void SelectedKinematicDecay::SetMissingQualityCriteria(const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction) {
  vtxSignPVRotSV_ = vtxSignPVRotSV;
  vtxSignPVRotPVRed_ = vtxSignPVRotPVRed;
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino                          
  //other channels need their own discriminators  
  //!!!
  a1Mass_ = a1Mass;
  energyTFraction_ = energyTFraction;
}


void SelectedKinematicDecay::SetInitalVertexProperties(reco::Vertex primaryVertexReFit,reco::Vertex primaryVertexReFitAndRotated,
						    std::vector<reco::TransientTrack> secVtxTracks, TransientVertex secVtx){
  initalPrimaryVertexReFit_=primaryVertexReFit;
  initalPrimaryVertexReFitAndRotated_=primaryVertexReFitAndRotated;
  for(unsigned int i=0;i<secVtxTracks.size();i++){
    initalSecVtxTracks_.push_back(secVtxTracks.at(i).track());
  }
  initalSecVtx_=secVtx;
}

void SelectedKinematicDecay::SetKFSecondaryVertex(reco::Vertex SecVtx){
  SecVtx_=SecVtx;
}

void SelectedKinematicDecay::SetInitalKinematics(TVector3 tauFlghtDir,std::vector<TLorentzVector> initalpions,TLorentzVector intial_a1_p4,double initThetaGJ, double ThetaMax){
  initalTauFlghtDir_=tauFlghtDir;
  initalpions_=initalpions;
  intial_a1_p4_=intial_a1_p4;
  initalThetaGJ_=initThetaGJ;
  initalThetaMax_=ThetaMax;
}

void SelectedKinematicDecay::SetInitalGuess(std::vector<TLorentzVector> &TauGuessLV,std::vector<TLorentzVector> &NuGuessLV){
  initalTauGuess_=TauGuessLV;
  initalNuGuess_=NuGuessLV;
}

TLorentzVector SelectedKinematicDecay::InitialTauGuess(unsigned int ambiguity){ 
  if(ambiguity<initalTauGuess_.size()) return initalTauGuess_.at(ambiguity);
  return TLorentzVector(0,0,0,0);
}

TLorentzVector SelectedKinematicDecay::InitalNeutrinoGuess(unsigned int ambiguity){
  if(ambiguity<initalNuGuess_.size())return initalNuGuess_.at(ambiguity);
  return TLorentzVector(0,0,0,0);
}

std::vector<TLorentzVector> SelectedKinematicDecay::InitalPions(){
  return  initalpions_;
}

TLorentzVector SelectedKinematicDecay::Inital_a1_p4(){
  /*  TLorentzVector intial_a1_p4(0,0,0,0);
  for(unsigned int i=0; i<initalpions_.size();i++){
    intial_a1_p4+=initalpions_.at(i);  
    }
    return intial_a1_p4;
  */
  return intial_a1_p4_;
}

TLorentzVector SelectedKinematicDecay::Tau(unsigned int ambiguity){
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="tau"){
      return iParticle->p4(); 
    }
  }
  return TLorentzVector(0,0,0,0); 
}

TLorentzVector SelectedKinematicDecay::Neutrino(unsigned int ambiguity){
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="neutrino"){
      return iParticle->p4();
    }
  }
  return TLorentzVector(0,0,0,0);
}

std::vector<TLorentzVector> SelectedKinematicDecay::Pions(unsigned int ambiguity){
  std::vector<TLorentzVector> pions;
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="pion"){
      pions.push_back(iParticle->p4());
    }
  }
  return pions;
}


TLorentzVector SelectedKinematicDecay::a1_p4(unsigned int ambiguity){
  TLorentzVector a1(0,0,0,0);
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="pion"){
      a1+=iParticle->p4();
    }
  }
  return a1;
}


const SelectedKinematicParticle * SelectedKinematicDecay::topParticle() const {
    return &(particles_.front());
}

void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle const * > & par) const {
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
  }
}

void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const {
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
      if ( std::abs(iter->charge()) == 1 ) {
	if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
      }
    }
}
void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const {
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if ( std::abs(iter->charge()) == 0 ) {
      if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
  }
}

void SelectedKinematicDecay::modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par) {
  for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if ( std::abs(iter->charge()) == 1 ) {
      if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
  }
}


const double SelectedKinematicDecay::chi2prob() const {
  ChiSquared chiSquared(chi2(), ndf());
  return chiSquared.probability();
}

const int SelectedKinematicDecay::sgnlConeTrkSize() const {
  if (PFTauRef().isAvailable()) {
    return PFTauRef()->signalPFChargedHadrCands().size();
  }
  printf("SelectedKinematicDecay::sgnlConeTrkSize:ERROR! reference to PFTau invalid. Return dummy value.\n");
  return 0;
}
