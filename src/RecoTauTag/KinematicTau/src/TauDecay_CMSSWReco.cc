#include "RecoTauTag/KinematicTau/interface/TauDecay_CMSSWReco.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"

#include <iomanip>
#include <cstdlib> 

TauDecay_CMSSWReco::TauDecay_CMSSWReco():
  TauDecay()
{

}

TauDecay_CMSSWReco::~TauDecay_CMSSWReco(){

} 
                            
bool TauDecay_CMSSWReco::AnalyzeTau(const reco::GenParticle *Tau,unsigned int &JAK_ID,unsigned int &TauBitMask){
  Reset();
  MotherIdx.clear();
  TauDecayProducts.clear();
  if(abs(Tau->pdgId())==PdtPdgMini::tau_minus){ // check that it is a tau
    unsigned int Tauidx=TauDecayProducts.size();
    TauDecayProducts.push_back(Tau);
    MotherIdx.push_back(Tauidx);
    for (unsigned int i=0; i< Tau->numberOfDaughters(); i++){
      const reco::Candidate *dau=Tau->daughter(i);
      Analyze(static_cast<const reco::GenParticle*>(dau),Tauidx);
    }
    ClassifyDecayMode(JAK_ID,TauBitMask);
    return true;
  }
  return false;
}
  
void TauDecay_CMSSWReco::Analyze(const reco::GenParticle *Particle,unsigned int midx){
  unsigned int pdgid=abs(Particle->pdgId());
  if(isTauFinalStateParticle(pdgid)){
    if(!isTauParticleCounter(pdgid)) std::cout << "TauDecay_CMSSWReco::Analyze WARNING: Unknown Final State Particle in Tau Decay... " << std::endl;
    TauDecayProducts.push_back(Particle);
    MotherIdx.push_back(midx);
    if(pdgid==PdtPdgMini::pi0 || PdtPdgMini::K_S0){// store information on pi0 || KS0 decay products even though they are final state particles 
      midx=MotherIdx.size()-1;
      for (unsigned int i=0; i< Particle->numberOfDaughters(); i++){
	const reco::Candidate *dau=Particle->daughter(i);
	AddFinalStateDecayInfo(static_cast<const reco::GenParticle*>(dau),midx);
      }
    }
    return;
  }
  if(Particle->status()==1 || isTauResonanceCounter(pdgid)){
    TauDecayProducts.push_back(Particle);
    MotherIdx.push_back(midx);
    midx=MotherIdx.size()-1;
  }
  for (unsigned int i=0; i< Particle->numberOfDaughters(); i++){
    const reco::Candidate *dau=Particle->daughter(i);
    Analyze(static_cast<const reco::GenParticle*>(dau),midx);
  }
}


void TauDecay_CMSSWReco::AddFinalStateDecayInfo(const reco::GenParticle *Particle,unsigned int midx){
  if(Particle->status()==1){
    TauDecayProducts.push_back(Particle);
    MotherIdx.push_back(midx);
    return;
  }
  for (unsigned int i=0; i< Particle->numberOfDaughters(); i++){
    const reco::Candidate *dau=Particle->daughter(i);
    AddFinalStateDecayInfo(static_cast<const reco::GenParticle*>(dau),midx);
  }
}


