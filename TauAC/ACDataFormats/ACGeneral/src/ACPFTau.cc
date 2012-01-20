#include "../interface/ACPFTau.h"

ClassImp(ACPFTau)

ACPFTau::ACPFTau() {
    p4_.SetXYZM(-1000,-1000,-1000,-1000);
    charge_ = -1000.0;
    discriminators_.clear();
    jetRef_ = ACJetRef();
    signalPFChargedHadrCands_ = -1;
    signalPFNeutrHadrCands_ = -1;
    signalPFGammaCands_ = -1;
    signalPiZeroCands_ = -1;
    isolationPFChargedHadrCands_ = -1;
    isolationPFNeutrHadrCands_ = -1;
    isolationPFGammaCands_ = -1;
    isolationPiZeroCands_ = -1;
}
ACPFTau::ACPFTau(const TLorentzVector& p4, int charge, const std::map<std::string, bool> & discriminators, int signalPFChargedHadrCands, int signalPFNeutrHadrCands, int signalPFGammaCands, int signalPiZeroCands, int isolationPFChargedHadrCands, int isolationPFNeutrHadrCands, int isolationPFGammaCands, int isolationPiZeroCands) {
    p4_ = p4;
    charge_ = charge;
    discriminators_ = discriminators;
    jetRef_ = ACJetRef();
    signalPFChargedHadrCands_ = signalPFChargedHadrCands;
    signalPFNeutrHadrCands_ = signalPFNeutrHadrCands;
    signalPFGammaCands_ = signalPFGammaCands;
    signalPiZeroCands_ = signalPiZeroCands;
    isolationPFChargedHadrCands_ = isolationPFChargedHadrCands;
    isolationPFNeutrHadrCands_ = isolationPFNeutrHadrCands;
    isolationPFGammaCands_ = isolationPFGammaCands;
    isolationPiZeroCands_ = isolationPiZeroCands;
}
const std::map<std::string, bool> & ACPFTau::discriminators() const { return discriminators_; }
const ACJetRef & ACPFTau::jetRef() const { return jetRef_;}
int ACPFTau::signalPFChargedHadrCands() const { return signalPFChargedHadrCands_; }
int ACPFTau::signalPFNeutrHadrCands() const { return signalPFNeutrHadrCands_; }
int ACPFTau::signalPFGammaCands() const { return signalPFGammaCands_; }
int ACPFTau::signalPiZeroCands() const { return signalPiZeroCands_; }
int ACPFTau::isolationPFChargedHadrCands() const { return isolationPFChargedHadrCands_; }
int ACPFTau::isolationPFNeutrHadrCands() const { return isolationPFNeutrHadrCands_; }
int ACPFTau::isolationPFGammaCands() const { return isolationPFGammaCands_; }
int ACPFTau::isolationPiZeroCands() const { return isolationPiZeroCands_; }
int ACPFTau::decayMode() const {
    // decay mode = 5*(number of charged hadrons-1) + number of pi0 candidates.
    if (signalPFChargedHadrCands()==0) return 0;
    int mode = 5 * (signalPFChargedHadrCands() - 1) + signalPiZeroCands();
    
    return mode;
}
