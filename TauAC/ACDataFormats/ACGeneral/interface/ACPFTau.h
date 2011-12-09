#ifndef ACPFTau_h
#define ACPFTau_h
/**
 Class to store jets.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <map>

#include "ACCandidate.h"

class ACPFTau : public ACCandidate {
public:
    virtual std::string classname() const { return "ACPFTau"; }
    ACPFTau();
    ACPFTau(const TLorentzVector& p4, int charge, const std::map<std::string, bool> & discriminators, int signalPFChargedHadrCands, int signalPFNeutrHadrCands);
    virtual ~ACPFTau() {};

    const std::map<std::string, bool> & discriminators() const;
    int signalPFChargedHadrCands() const;
    int signalPFNeutrHadrCands() const;

    ClassDef(ACPFTau,1);

protected:
    /// results of the original PFTau discriminators
    std::map<std::string, bool> discriminators_;
    /// size of tracks in the signal cone of the PFTau candidate
    int signalPFChargedHadrCands_;
    /// size of neutral candidates in the signal cone of the PFTau candidate
    int signalPFNeutrHadrCands_;
};

typedef ACPFTau* pACPFTau;

#endif
