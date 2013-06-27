#include "../interface/ACMuon.h"

ClassImp(ACMuon)

ACMuon::ACMuon() {
    p4_.SetXYZM(-1000,-1000,-1000,-1000);
    vertex_.SetXYZ(-1000, -1000, -1000);
    charge_ = -1000.0;
    name_ = "";
    
    trackRef_ = ACTrackRef ();
    genRef_ = ACGenParticleRef();
    isGlobalMuonPromptTight_ = false;
    isRecommendedMuon_ = false;
    sumChargedHadronPt_ = -1.;
    sumChargedParticlePt_ = -1.;
    sumNeutralHadronEt_ = -1.;
    sumPhotonEt_ = -1.;
    sumPUPt_ = -1.;
}
ACMuon::ACMuon(const TLorentzVector& p4, int charge) {
    p4_ = p4;
    charge_ = charge;
    vertex_.SetXYZ(-1000, -1000, -1000);
    name_ = "muon";
    
    trackRef_ = ACTrackRef ();
    genRef_ = ACGenParticleRef();
    isGlobalMuonPromptTight_ = false;
    isRecommendedMuon_ = false;
    sumChargedHadronPt_ = -1.;
    sumChargedParticlePt_ = -1.;
    sumNeutralHadronEt_ = -1.;
    sumPhotonEt_ = -1.;
    sumPUPt_ = -1.;
}

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
double ACMuon::vx() const { return vertex().X(); }
double ACMuon::vy() const { return vertex().Y(); }
double ACMuon::vz() const { return vertex().Z(); }
double ACMuon::isolation() const {
    double value = sumChargedHadronPt() + sumNeutralHadronEt() + sumPhotonEt();
    if (pt() > 0. ) return value/pt();
    else return INFINITY;
}
double ACMuon::isolationDBeta() const {
    double value = sumChargedHadronPt() + std::max(0., sumNeutralHadronEt() + sumPhotonEt() - .5 * sumPUPt());
    if (pt() > 0. ) return value/pt();
    else return INFINITY;
}
bool ACMuon::isIsolatedLoose() const {
    return isolation() <= 0.2;
}
bool ACMuon::isIsolatedTight() const {
    return isolation() <= 0.12;
}
