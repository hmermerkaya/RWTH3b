#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

SelectedKinematicDecay::SelectedKinematicDecay() {
  SetInitialProperties(Undefined,reco::PFTauRef(),std::vector<reco::TrackRef>(),reco::Vertex(),"",0);
  SetKinematicFitProperties(NAmbiguity,SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0);
  SetKinematicFitStatus(NAmbiguity,std::map<std::string, bool>());
  SetQualityCriteria(NAmbiguity,-1,-1,-1,-1);
  SetInitialVertexProperties(reco::Vertex(),reco::Vertex(),std::vector<reco::TransientTrack>(),TransientVertex());  
  SetInitialKinematics(TVector3(),std::vector<TLorentzVector>(),TLorentzVector(),TVector3(),0.0,0.0);
}

SelectedKinematicDecay::SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRefOrig, std::vector<reco::TrackRef> &TrackTriplet,
					       const reco::Vertex &primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx){
  SetInitialProperties(tauDecayMode,tauRefOrig,TrackTriplet,primaryVertex,primVtxReFitTag,nTauPerVtx);
  SetKinematicFitProperties(NAmbiguity,SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0);
  SetKinematicFitStatus(NAmbiguity,std::map<std::string, bool>());
  SetQualityCriteria(NAmbiguity,-1,-1,-1,-1);
  SetInitialVertexProperties(reco::Vertex(),reco::Vertex(),std::vector<reco::TransientTrack>(),TransientVertex());
  SetInitialKinematics(TVector3(),std::vector<TLorentzVector>(),TLorentzVector(),TVector3(),0.0,0.0);
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

  initialTrackTriplet_=TrackTriplet;
  initialPrimVtx_=primaryVertex;
}

void SelectedKinematicDecay::SetKinematicFitStatus(unsigned int ambiguity, const std::map<std::string, bool> discriminators){
  if(discriminators_.size()!=NAmbiguity){discriminators_.clear();discriminators_.resize(NAmbiguity);}
  if(ambiguity<NAmbiguity) discriminators_.at(ambiguity) = discriminators;
}

void SelectedKinematicDecay::SetKinematicFitProperties(unsigned int ambiguity, const SelectedKinematicParticleCollection particles,
						       const int iterations, const int maxiterations, const float csum,
						       const float mincsum, const int constraints, const int ndf,
						       const float chi2){


  if(iterations_.size()!=NAmbiguity){iterations_.clear();iterations_.resize(NAmbiguity,-1);}
  if(maxiterations_.size()!=NAmbiguity){maxiterations_.clear();maxiterations_.resize(NAmbiguity,-1);}
  if(csum_.size()!=NAmbiguity){csum_.clear();csum_.resize(NAmbiguity,-1);}
  if(mincsum_.size()!=NAmbiguity){mincsum_.clear();mincsum_.resize(NAmbiguity,-1);}
  if(constraints_.size()!=NAmbiguity){constraints_.clear();constraints_.resize(NAmbiguity,-1);}
  if(ndf_.size()!=NAmbiguity){ndf_.clear();ndf_.resize(NAmbiguity,0);}
  if(chi2_.size()!=NAmbiguity){chi2_.clear();chi2_.resize(NAmbiguity,0);}

  if(ambiguity<NAmbiguity){
    //clean particles with the same ambiguity
    for(unsigned int i=0;i<particles_.size();i++){
      if(particles_.at(i).ambiguity()==ambiguity){
	particles_.erase(particles_.begin()+i);
	i--;
      }
    }
    //and new particles 
    for(unsigned int i=0;i<particles.size();i++){
      particles_.push_back(particles.at(i));
    }
    iterations_.at(ambiguity) = iterations;
    maxiterations_.at(ambiguity) = maxiterations;
    csum_.at(ambiguity) = csum;
    mincsum_.at(ambiguity) = mincsum;
    constraints_.at(ambiguity) = constraints;
    ndf_.at(ambiguity) = ndf;
    chi2_.at(ambiguity) = chi2;
  }
}


void SelectedKinematicDecay::SetQualityCriteria(unsigned int ambiguity,const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction) {
  if(vtxSignPVRotSV_.size()!=NAmbiguity){vtxSignPVRotSV_.clear();vtxSignPVRotSV_.resize(NAmbiguity,-1);}
  if(vtxSignPVRotPVRed_.size()!=NAmbiguity){vtxSignPVRotPVRed_.clear();vtxSignPVRotPVRed_.resize(NAmbiguity,-1);}
  if(a1Mass_.size()!=NAmbiguity){a1Mass_.clear();a1Mass_.resize(NAmbiguity,-1);}
  if(energyTFraction_.size()!=NAmbiguity){energyTFraction_.clear();energyTFraction_.resize(NAmbiguity,-1);}
  if(ambiguity<NAmbiguity){
    vtxSignPVRotSV_.at(ambiguity) = vtxSignPVRotSV;
    vtxSignPVRotPVRed_.at(ambiguity) = vtxSignPVRotPVRed;
    //WARNING!!!
    //from now one we assume a tau decay into three pions and neutrino                          
    //other channels need their own discriminators  
    //!!!
    a1Mass_.at(ambiguity) = a1Mass;
    energyTFraction_.at(ambiguity) = energyTFraction;
  }
}


void SelectedKinematicDecay::SetInitialVertexProperties(reco::Vertex primaryVertexReFit,reco::Vertex primaryVertexReFitAndRotated,
						    std::vector<reco::TransientTrack> secVtxTracks, TransientVertex secVtx){
  initialPrimaryVertexReFit_=primaryVertexReFit;
  initialPrimaryVertexReFitAndRotated_=primaryVertexReFitAndRotated;
  for(unsigned int i=0;i<secVtxTracks.size();i++){
    initialSecVtxTracks_.push_back(secVtxTracks.at(i).track());
  }
  initialSecVtx_=secVtx;
}

void SelectedKinematicDecay::SetKFSecondaryVertex(unsigned int ambiguity,reco::Vertex SecVtx){
  if(SecVtx_.size()!=NAmbiguity){SecVtx_.clear();SecVtx_.resize(NAmbiguity);}
  if(ambiguity<NAmbiguity)SecVtx_.at(ambiguity)=SecVtx;
}

void SelectedKinematicDecay::SetInitialKinematics(TVector3 tauFlghtDirNoCorr,std::vector<TLorentzVector> initialpions,
						 TLorentzVector initial_a1_p4,TVector3 tauFlghtDir,double initThetaGJ, double ThetaMax){
  initialTauFlghtDirNoCorr_=tauFlghtDirNoCorr;
  initialTauFlghtDir_=tauFlghtDir;
  initialpions_=initialpions;
  initial_a1_p4_=initial_a1_p4;
  initialThetaGJ_=initThetaGJ;
  initialThetaMax_=ThetaMax;
}

void SelectedKinematicDecay::SetInitialGuess(unsigned int ambiguity,TLorentzVector &TauGuessLV,TLorentzVector &NuGuessLV,TVector3 &TauFlghtDirGuess){
  if(initialTauGuess_.size()!=NAmbiguity){initialTauGuess_.clear();initialTauGuess_.resize(NAmbiguity);}
  if(initialNuGuess_.size()!=NAmbiguity){initialNuGuess_.clear();initialNuGuess_.resize(NAmbiguity);}
  if(TauFlghtDirGuess_.size()!=NAmbiguity){TauFlghtDirGuess_.clear();TauFlghtDirGuess_.resize(NAmbiguity);}
  if(ambiguity<NAmbiguity){
    initialTauGuess_.at(ambiguity)=TauGuessLV;
    initialNuGuess_.at(ambiguity)=NuGuessLV;
    TauFlghtDirGuess_.at(ambiguity)=TauFlghtDirGuess;
  }
}

TLorentzVector SelectedKinematicDecay::InitialTauGuess(unsigned int ambiguity){ 
  if(ambiguity<initialTauGuess_.size()) return initialTauGuess_.at(ambiguity);
  return TLorentzVector(0,0,0,0);
}

TLorentzVector SelectedKinematicDecay::InitialNeutrinoGuess(unsigned int ambiguity){
  if(ambiguity<initialNuGuess_.size())return initialNuGuess_.at(ambiguity);
  return TLorentzVector(0,0,0,0);
}

TVector3 SelectedKinematicDecay::InitialTauFlghtDirGuess(unsigned int ambiguity){
  if(ambiguity<TauFlghtDirGuess_.size())return TauFlghtDirGuess_.at(ambiguity);
  return TVector3(0,0,0);
}

std::vector<TLorentzVector> SelectedKinematicDecay::InitialPions(){
  return  initialpions_;
}

TLorentzVector SelectedKinematicDecay::Initial_a1_p4(){
  return initial_a1_p4_;
}

TLorentzVector SelectedKinematicDecay::Tau(unsigned int ambiguity){
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="tau"  && iParticle->ambiguity()==ambiguity){
      return iParticle->p4(); 
    }
  }
  return TLorentzVector(0,0,0,0); 
}

TLorentzVector SelectedKinematicDecay::Neutrino(unsigned int ambiguity){
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="neutrino"  && iParticle->ambiguity()==ambiguity){
      return iParticle->p4();
    }
  }
  return TLorentzVector(0,0,0,0);
}

std::vector<TLorentzVector> SelectedKinematicDecay::Pions(unsigned int ambiguity){
  std::vector<TLorentzVector> pions;
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="pion"  && iParticle->ambiguity()==ambiguity){
      pions.push_back(iParticle->p4());
    }
  }
  return pions;
}


TLorentzVector SelectedKinematicDecay::a1_p4(unsigned int ambiguity){
  TLorentzVector a1(0,0,0,0);
  for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = particles_.begin(); iParticle != particles_.end(); ++iParticle){
    if(iParticle->name()=="pion" && iParticle->ambiguity()==ambiguity){
      a1+=iParticle->p4();
    }
  }
  return a1;
}


const SelectedKinematicParticle * SelectedKinematicDecay::topParticle(unsigned int ambiguity) const {
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if(iter->name()=="tau" && iter->ambiguity()==ambiguity) return &(*iter);
  }
  return &(*particles_.begin());
}

void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity) const {
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if(iter->name()!="tau" && iter->ambiguity()==ambiguity) par.push_back(&(*iter));//skip mother
  }
}

void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity) const {
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
      if ( std::abs(iter->charge()) == 1 && iter->ambiguity()==ambiguity) {
	if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
      }
    }
}
void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity) const {
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if ( std::abs(iter->charge()) == 0 && iter->ambiguity()==ambiguity ) {
      if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
  }
}

void SelectedKinematicDecay::modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par,unsigned int ambiguity) {
  for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if ( std::abs(iter->charge()) == 1  ) {
      if(iter != particles_.begin() && iter->name()!="tau"  && iter->name()!="a1") par.push_back(&(*iter));//skip mother
    }
  }
}


const double SelectedKinematicDecay::chi2prob(unsigned int ambiguity) const {
  ChiSquared chiSquared(chi2(ambiguity), ndf(ambiguity));
  return chiSquared.probability();
}

const int SelectedKinematicDecay::sgnlConeTrkSize() const {
  if (PFTauRef().isAvailable()) {
    return PFTauRef()->signalPFChargedHadrCands().size();
  }
  printf("SelectedKinematicDecay::sgnlConeTrkSize:ERROR! reference to PFTau invalid. Return dummy value.\n");
  return 0;
}

reco::Vertex SelectedKinematicDecay::SecondaryVertex(unsigned int ambiguity){
  if(ambiguity<NAmbiguity) return SecVtx_.at(ambiguity);
  return reco::Vertex();
}

// KF variables                                               
const std::map<std::string, bool> SelectedKinematicDecay::discriminators(unsigned int ambiguity)const{ 
  if(ambiguity<NAmbiguity) return discriminators_.at(ambiguity);
  return std::map<std::string, bool>(); 
}

const int SelectedKinematicDecay::iterations(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return iterations_.at(ambiguity);
  return -999;
}

const int SelectedKinematicDecay::maxiterations(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return maxiterations_.at(ambiguity);
  return -999;
}

const float  SelectedKinematicDecay::chi2(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return chi2_.at(ambiguity);
  return -999;
}

const float SelectedKinematicDecay::constraints(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return constraints_.at(ambiguity);
  return -999;
}

const float SelectedKinematicDecay::ndf(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return ndf_.at(ambiguity);
  return -999;
}

const float SelectedKinematicDecay::csum(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return csum_.at(ambiguity);
  return -999;
}

const float SelectedKinematicDecay::mincsum(unsigned int ambiguity)const{
  if(ambiguity<NAmbiguity) return mincsum_.at(ambiguity);
  return -999;
}

const double SelectedKinematicDecay::vtxSignPVRotSV(unsigned int ambiguity)const{ 
  if(ambiguity<NAmbiguity) return vtxSignPVRotSV_.at(ambiguity); 
  return -999;
}  

const double SelectedKinematicDecay::vtxSignPVRotPVRed(unsigned int ambiguity)const{ 
  if(ambiguity<NAmbiguity) return vtxSignPVRotPVRed_.at(ambiguity);
  return -999;
}

const double SelectedKinematicDecay::a1Mass(unsigned int ambiguity)const{ 
  if(ambiguity<NAmbiguity) return a1Mass_.at(ambiguity);
  return -999;
}             
          
const double SelectedKinematicDecay::energyTFraction(unsigned int ambiguity)const{ 
  if(ambiguity<NAmbiguity) return energyTFraction_.at(ambiguity); 
  return -999;
}    


const SelectedKinematicParticleCollection SelectedKinematicDecay::particles(unsigned int ambiguity) const{
  SelectedKinematicParticleCollection par;
  for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
    if(iter->ambiguity()==ambiguity) par.push_back((*iter));
  }
  return par;
}
