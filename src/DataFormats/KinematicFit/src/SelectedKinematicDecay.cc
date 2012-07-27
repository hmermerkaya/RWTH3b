#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

SelectedKinematicDecay::SelectedKinematicDecay() {
  SetInitialProperties(Undefined,reco::PFTauRef(),std::vector<reco::TrackRef>(),reco::Vertex(),"",0);
  SetKinematicFitProperties(SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0,std::map<std::string, bool>());
  setMissingQualityCriteria(-1,-1,-1,-1);
  SetPrimaryVertexReFit(reco::Vertex());
  SetPrimaryVertexReFitAndRotated(reco::Vertex());
  SetSecondaryVertex(std::vector<reco::TransientTrack>(),TransientVertex());
}


SelectedKinematicDecay::SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRef, std::vector<reco::TrackRef> &TrackTriplet,
					       const reco::Vertex &primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx){
  SetInitialProperties(tauDecayMode,tauRef,TrackTriplet,primaryVertex,primVtxReFitTag,nTauPerVtx);
  SetKinematicFitProperties(SelectedKinematicParticleCollection(),-1,-1,-1.0,-1.0,-1,0,0.0,std::map<std::string, bool>());
  setMissingQualityCriteria(-1,-1,-1,-1);
  SetPrimaryVertexReFit(reco::Vertex());
  SetPrimaryVertexReFitAndRotated(reco::Vertex());
  SetSecondaryVertex(std::vector<reco::TransientTrack>(),TransientVertex());
}

SelectedKinematicDecay::SelectedKinematicDecay(const SelectedKinematicParticleCollection & particles,
						 const int iterations, const int maxiterations, const float csum,
						 const float mincsum, const int constraints, const int ndf, const float chi2,
						 const reco::PFTauRef & tauRef, const std::map<std::string, bool> & discriminators){
  SetInitialProperties(Undefined,tauRef,std::vector<reco::TrackRef>(),reco::Vertex()," ",0);
  SetKinematicFitProperties(particles,iterations,maxiterations,csum,mincsum,constraints,ndf,chi2,discriminators);
  setMissingQualityCriteria(-1,-1,-1,-1);
  SetPrimaryVertexReFit(reco::Vertex());
  SetPrimaryVertexReFitAndRotated(reco::Vertex()); 
  SetSecondaryVertex(std::vector<reco::TransientTrack>(),TransientVertex());
}

SelectedKinematicDecay::~SelectedKinematicDecay(){

}

void SelectedKinematicDecay::SetInitialProperties(unsigned int tauDecayMode,reco::PFTauRef tauRef,std::vector<reco::TrackRef> TrackTriplet, const reco::Vertex primaryVertex,
						  std::string primVtxReFitTag, unsigned int nTauPerVtx){
  tauDecayMode_=tauDecayMode;
  PFTauRef_=tauRef;
  TrackTriplet_=TrackTriplet;
  primVtx_=primaryVertex;
  primVtxReFitTag_=primVtxReFitTag;
  nTauPerVtx_=nTauPerVtx;
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

void SelectedKinematicDecay::setMissingQualityCriteria(const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction) {
  vtxSignPVRotSV_ = vtxSignPVRotSV;
  vtxSignPVRotPVRed_ = vtxSignPVRotPVRed;
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino                          
  //other channels need their own discriminators  
  //!!!
  a1Mass_ = a1Mass;
  energyTFraction_ = energyTFraction;
}


void SelectedKinematicDecay::SetPrimaryVertexReFit(reco::Vertex primaryVertexReFit){
  primaryVertexReFit_=primaryVertexReFit;
}

void SelectedKinematicDecay::SetPrimaryVertexReFitAndRotated(reco::Vertex primaryVertexReFitAndRotated){
  primaryVertexReFitAndRotated_=primaryVertexReFitAndRotated;
}

void SelectedKinematicDecay::SetSecondaryVertex(std::vector<reco::TransientTrack> secVtxTracks, TransientVertex secVtx){
  for(unsigned int i=0;i<secVtxTracks.size();i++){
    secVtxTracks_.push_back(secVtxTracks.at(i).track());
  }
  secVtx_=secVtx;
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

const int SelectedKinematicDecay::iterations() const {
    return iterations_;
}
const int SelectedKinematicDecay::maxiterations() const {
    return maxiterations_;
}
const float SelectedKinematicDecay::chi2() const {
    return chi2_;
}
const float SelectedKinematicDecay::constraints() const {
    return constraints_;
}
const float SelectedKinematicDecay::ndf() const {
    return ndf_;
}
const float SelectedKinematicDecay::csum() const {
    return csum_;
}
const float SelectedKinematicDecay::mincsum() const {
    return mincsum_;
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
