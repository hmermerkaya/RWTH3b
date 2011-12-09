#ifndef ACFittedDecay_H
#define ACFittedDecay_H
/**
 Derived class for tau decays created by a kinematic fit.

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */

//#include <map>

#include "ACDecay.h"
#include "ACFitParticle.h"

class ACFittedDecay : public ACDecayBase<ACFitParticle> {
public:
    ACFittedDecay();
    ACFittedDecay(std::vector<ACFitParticleRef> particles, std::string name, int iterations, float csum, int maxiterations, float mincsum, float chi2, ACPFTauRef PFTauRef);
    virtual ~ACFittedDecay() {};

    int iterations() const;
    float csum() const;
    int maxiterations() const;
    float mincsum() const;
    float chi2() const;
    const ACPFTauRef & PFTauRef() const;

    /// return the initial primary vertex of the tau decay
    // ACVertex initialPrimaryVertex() const;

    ClassDef(ACFittedDecay,1);

protected:
    /// iterations of the kinematic fit
    int iterations_;
    /// sum of constraints of the kinematic fit
    float csum_;
    /// maximum of allowed iterations (static parameter of the kinematic fit)
    int maxiterations_;
    /// minimal sum of constraints (static parameter of the kinematic fit)
    float mincsum_;
    /// chi2 of the kinematic fit (valid ndf only available in derived class with a decay mode assumption)
    float chi2_;
    /// reference to the initial PFTau from which this decay was created
    ACPFTauRef PFTauRef_;
};

typedef ACFittedDecay* pACFittedDecay;

#endif
