#ifndef ACPFTau_h
#define ACPFTau_h
/**
 Class to store offline PFTaus.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <map>

#include "ACCandidate.h"
#include "ACRef.h"

class ACPFTau : public ACCandidate {
public:
    virtual std::string classname() const { return "ACPFTau"; }
    ACPFTau();
    ACPFTau(const TLorentzVector& p4, int charge, const std::map<std::string, bool> & discriminators, int signalPFChargedHadrCands, int signalPFNeutrHadrCands, int signalPFGammaCands, int signalPiZeroCands, int isolationPFChargedHadrCands, int isolationPFNeutrHadrCands, int isolationPFGammaCands, int isolationPiZeroCands);
    virtual ~ACPFTau() {};

    const std::map<std::string, bool> & discriminators() const;
    const ACJetRef & jetRef() const;
    int signalPFChargedHadrCands() const;
    int signalPFNeutrHadrCands() const;
    int signalPFGammaCands() const;
    int signalPiZeroCands() const;
    int isolationPFChargedHadrCands() const;
    int isolationPFNeutrHadrCands() const;
    int isolationPFGammaCands() const;
    int isolationPiZeroCands() const;

    /// returns the decay mode according to the TauPOG definition: mode = 5 * (chargedDaughters - 1) + number of pi0's
    int decayMode() const;
    ClassDef(ACPFTau,1);

protected:
    /// results of the original PFTau discriminators
    std::map<std::string, bool> discriminators_;
    /// reference to initial jet (FIXME: not filled yet)
    ACJetRef jetRef_;

    /// size of tracks in the signal cone of the PFTau candidate
    int signalPFChargedHadrCands_;
    /// size of neutral candidates in the signal cone of the PFTau candidate
    int signalPFNeutrHadrCands_;
    /// size of gamma candidates in the signal cone of the PFTau candidate
    int signalPFGammaCands_;
    /// Association of gamma candidates into pi0's in the signal cone of the PFTau candidate
    int signalPiZeroCands_;

    /// size of tracks in the isolation cone of the PFTau candidate
    int isolationPFChargedHadrCands_;
    /// size of neutral candidates in the isolation cone of the PFTau candidate
    int isolationPFNeutrHadrCands_;
    /// size of gamma candidates in the isolation cone of the PFTau candidate
    int isolationPFGammaCands_;
    /// Association of gamma candidates into pi0's in the isolation cone of the PFTau candidate
    int isolationPiZeroCands_;
};

typedef ACPFTau* pACPFTau;

#endif
