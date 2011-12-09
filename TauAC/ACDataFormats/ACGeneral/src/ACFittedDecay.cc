#include "../interface/ACFittedDecay.h"

ClassImp(ACFittedDecay)

ACFittedDecay::ACFittedDecay():
ACDecayBase<ACFitParticle>() {
    iterations_ = -1;
    csum_ = -1.;
    maxiterations_ = -1;
    mincsum_ = -1.;
    PFTauRef_ = ACPFTauRef();
}
ACFittedDecay::ACFittedDecay(std::vector<ACRef<ACFitParticle> > particles, std::string name, int iterations, float csum, int maxiterations, float mincsum, float chi2, ACPFTauRef PFTauRef):
ACDecayBase<ACFitParticle>(particles, name) {
    iterations_ = iterations;
    csum_ = csum;
    maxiterations_ = maxiterations;
    mincsum_ = mincsum;
    chi2_ = chi2;
    PFTauRef_ = PFTauRef;
}

int ACFittedDecay::iterations() const { return iterations_; }
float ACFittedDecay::csum() const { return csum_; }
int ACFittedDecay::maxiterations() const { return maxiterations_; }
float ACFittedDecay::mincsum() const { return mincsum_; }
float ACFittedDecay::chi2() const { return chi2_; }
const ACPFTauRef & ACFittedDecay::PFTauRef() const { return PFTauRef_; }
