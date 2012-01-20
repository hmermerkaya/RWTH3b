#ifndef ACFittedThreeProngDecay_H
#define ACFittedThreeProngDecay_H
/**
 Derived class for 3-prong tau decays created by a kinematic fit.

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */

#include "ACFittedDecay.h"

class ACFittedThreeProngDecay : public ACFittedDecay {
public:
    virtual std::string classname() const { return "ACFittedThreeProngDecay"; }
    ACFittedThreeProngDecay();
    ACFittedThreeProngDecay(std::vector<ACRef<ACFitParticle> > particles, std::string name, int iterations, float csum, int maxiterations, float mincsum, float chi2, ACPFTauRef PFTauRef, int constraints);
    virtual ~ACFittedThreeProngDecay() {};

    /**
     This is the real ndf derived from the KinematicFit (and not what is called ndf in the KinematicFit).
     It is the number of constraints used in a particular decay minus the number of initial free parameters (under a specific decay assumption).
     */
    int ndf() const;
    /// normalized chi2 of the kinematic fit
    float normalizedChi2() const;
    /// chi2 probability of the kinematic fit
    float chi2Probability() const;

    const ACFitParticle * tau() const;
    const ACFitParticle * prong0() const;
    const ACFitParticle * prong1() const;
    const ACFitParticle * prong2() const;
    const ACFitParticle * nu() const;

    /**
     returns summed lorentz vector of all visible daughters. this is the a1 for this specific decay assumption. per default the fitted state is regarded.
     */
    TLorentzVector a1Lorentz(bool prefit = false) const;
    /// returns the GJ angle. per default the fitted state is regarded.
    double gjAngle(bool prefit = false) const;
    /// returns the maximal allowed GJ angle. per default the fitted state is regarded.
    double maxGJAngle(bool prefit = false) const;

    ClassDef(ACFittedThreeProngDecay,1);

protected:
    /// constraints of the kinematic fit (called ndf in the KinematicFit package)
    int constraints_;

private:
    /// integrity check for correct amount of particles
    bool testDecay() const;
};

typedef ACFittedThreeProngDecay* pACFittedThreeProngDecay;

#endif
