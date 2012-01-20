#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"



SelectedKinematicDecay::SelectedKinematicDecay() {

}
SelectedKinematicDecay::SelectedKinematicDecay(SelectedKinematicParticleCollection particles)
{
    particles_ = particles;
	signalPFChargedHadrCands_ = -1;
	signalPFNeutrHadrCands_ = -1;
	fraction_ = -1;	       
        RefitMass_ = -1;	       
        Chi2_ = -1;		       
        ModPV_PV_significance_ = -1;
	SV_PV_significance_ = -1;   
	discriminators_.clear();
	PFTauRef_.clear();
}
SelectedKinematicDecay::SelectedKinematicDecay(SelectedKinematicParticleCollection particles, const int & signalPFChargedHadrCands, const int & signalPFNeutrHadrCands, 
					       const double & fraction, const double & RefitMass, const double & Chi2, const double & ModPV_PV_significance,const double & SV_PV_significance, const std::map<std::string, bool> & discriminators)
{
    particles_ = particles;
	signalPFChargedHadrCands_ = signalPFChargedHadrCands;
	signalPFNeutrHadrCands_   = signalPFNeutrHadrCands;
	discriminators_ = discriminators;
	fraction_ = fraction;	       
        RefitMass_ = RefitMass;	       
        Chi2_ = Chi2;		       
        ModPV_PV_significance_ = ModPV_PV_significance;
	SV_PV_significance_ = SV_PV_significance;   

	PFTauRef_.clear();
}

const SelectedKinematicParticle* SelectedKinematicDecay::topParticle() const
{
    return &(particles_.front());
}
void SelectedKinematicDecay::particles(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		par.push_back(&(*iter));
    }
}
void SelectedKinematicDecay::daughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
		if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
    }
}
void SelectedKinematicDecay::chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
void SelectedKinematicDecay::neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const
{
    for ( SelectedKinematicParticleCollection::const_iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 0 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}
int SelectedKinematicDecay::signalPFChargedHadrCands() const {
	return signalPFChargedHadrCands_;
}
int SelectedKinematicDecay::signalPFNeutrHadrCands() const {
	return signalPFNeutrHadrCands_;
}

double SelectedKinematicDecay::TauEnergyFraction() const {
  return fraction_;
}
double SelectedKinematicDecay::RefitVisibleMass() const {
  return RefitMass_;
}

double SelectedKinematicDecay::Chi2() const {
return Chi2_;
}

double SelectedKinematicDecay::ModPV_PV_significance() const {
  return ModPV_PV_significance_;
}

double SelectedKinematicDecay::PV_SV_significance() const {
  return SV_PV_significance_;
}







std::map<std::string, bool> SelectedKinematicDecay::discriminators() const
{
	return discriminators_;
}
void SelectedKinematicDecay::modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par)
{
    for ( SelectedKinematicParticleCollection::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        if ( std::abs(iter->charge()) == 1 ) {
			if(iter != particles_.begin()) par.push_back(&(*iter));//skip mother
        }
    }
}

void SelectedKinematicDecay::setPFTauRef(const std::vector<reco::PFTauRef> & value){
	PFTauRef_ = value;
}
void SelectedKinematicDecay::addPFTauRef(const reco::PFTauRef & value){
	PFTauRef_.push_back(value);
}
std::vector<reco::PFTauRef> SelectedKinematicDecay::PFTauRef() const{
	return PFTauRef_;
}
