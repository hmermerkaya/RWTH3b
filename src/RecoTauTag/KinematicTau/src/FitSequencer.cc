#include "RecoTauTag/KinematicTau/interface/FitSequencer.h"

FitSequencer::FitSequencer(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder,edm::Handle<reco::GenParticleCollection> &GenPart):
  transientTrackBuilder_(transTrackBuilder),
  GenPart_(GenPart),
  status(0)
{

}
FitSequencer::FitSequencer(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const edm::ParameterSet& cfg,edm::Handle<reco::GenParticleCollection> &GenPart):
  transientTrackBuilder_(transTrackBuilder),
  GenPart_(GenPart),
  status(0)
{

}


void FitSequencer::GetSelectedKinematicParticleList(unsigned int &ambiguity,SelectedKinematicParticleCollection &refitTauDecay){
  refitTauDecay.clear();
  reco::RecoChargedCandidateRef emptyCandRef;// temporary dummy variable
  refitTauDecay.push_back(SelectedKinematicParticle(mother(),status,ambiguity, emptyCandRef));
  std::vector<LorentzVectorParticle> d=daugthers(); 
  for(unsigned int i=0;i<d.size();i++){
    refitTauDecay.push_back(SelectedKinematicParticle(d.at(i),status,ambiguity,emptyCandRef));
  }
}

void FitSequencer::StoreResults(float chi2,float ndf,float csum,float niter,float nconst,std::vector<LorentzVectorParticle> fitdaughters,LorentzVectorParticle fitmother){
  Chi2_.push_back(chi2);
  NDF_.push_back(ndf);
  cSum_.push_back(csum);
  Niter_.push_back(niter);
  NConstraints_.push_back(nconst);
  FitMother_.push_back(fitmother);
  FitDaughters_.push_back(fitdaughters);
}

LorentzVectorParticle FitSequencer::mother(int i){
  if(FitMother_.size()==0) return LorentzVectorParticle();
  if(i<0 || i>=(int)FitMother_.size()) i=FitMother_.size()-1;
  return FitMother_.at(i);
}

std::vector<LorentzVectorParticle> FitSequencer::daugthers(int i){
  if(FitDaughters_.size()==0) return std::vector<LorentzVectorParticle>();
  if(i<0){
    std::vector<LorentzVectorParticle> d;
    for(unsigned int i=0;i<FitDaughters_.size();i++){
      for(unsigned int j=0;j<FitDaughters_.at(i).size();j++)d.push_back(FitDaughters_.at(i).at(j));
    }
    return d;
  }
  if(i>=(int)FitDaughters_.size()) i=FitDaughters_.size()-1;
  return FitDaughters_.at(i);
}

std::vector<LorentzVectorParticle> FitSequencer::chargedDaughters(int i){
  std::vector<LorentzVectorParticle> cd;
  std::vector<LorentzVectorParticle> d=daugthers(i);
  for(unsigned int i=0;i<d.size();i++){if(fabs(d.at(i).Charge()>0.1))cd.push_back(d.at(i));}
  return cd;
}

std::vector<LorentzVectorParticle> FitSequencer::neutralDaughters(int i){
  std::vector<LorentzVectorParticle> nd;
  std::vector<LorentzVectorParticle> d=daugthers(i);
  for(unsigned int i=0;i<d.size();i++){if(fabs(d.at(i).Charge()<0.1))nd.push_back(d.at(i));}
  return nd;
}

reco::PFTau FitSequencer::getPFTau(){
  TLorentzVector a1(A1.LV());
  TLorentzVector tau(Tau.LV());
  math::XYZTLorentzVector ptau(tau.Px(),tau.Py(),tau.Pz(),tau.E());
  math::XYZTLorentzVector pa1(a1.Px(),a1.Py(),a1.Pz(),a1.E());
  reco::PFTau refitPFTau(int(Tau.Charge()), pa1, PV_.position());
  refitPFTau.setalternatLorentzVect(ptau);
  return refitPFTau;
}
