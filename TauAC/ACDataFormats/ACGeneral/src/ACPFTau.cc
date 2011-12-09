#include "../interface/ACPFTau.h"

ClassImp(ACPFTau)

ACPFTau::ACPFTau() {
    p4_.SetXYZM(-1000,-1000,-1000,-1000);
    charge_ = -1000.0;
    discriminators_.clear();
    signalPFChargedHadrCands_ = -1;
    signalPFNeutrHadrCands_ = -1;
}
ACPFTau::ACPFTau(const TLorentzVector& p4, int charge, const std::map<std::string, bool> & discriminators, int signalPFChargedHadrCands, int signalPFNeutrHadrCands) {
    p4_ = p4;
    charge_ = charge;
    discriminators_ = discriminators;
    signalPFChargedHadrCands_ = signalPFChargedHadrCands;
    signalPFNeutrHadrCands_ = signalPFNeutrHadrCands;
}
int ACPFTau::signalPFChargedHadrCands() const { return signalPFChargedHadrCands_; }
int ACPFTau::signalPFNeutrHadrCands() const { return signalPFNeutrHadrCands_; }
const std::map<std::string, bool> & ACPFTau::discriminators() const { return discriminators_; }
