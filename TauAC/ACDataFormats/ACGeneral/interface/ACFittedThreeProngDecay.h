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
    ACFittedThreeProngDecay(std::vector<ACFitParticleRef> particles, std::string name, int iterations, float csum, int maxiterations, float mincsum, float chi2, ACPFTauRef PFTauRef, int constraints, int ndf);
    virtual ~ACFittedThreeProngDecay() {};

    /**
     This is the real ndf derived from the KinematicFit (and not what is called ndf in the RecoVertex/KinematicFit).
     It is the number of constraints used in a particular decay minus the number of initial free parameters (under a specific decay assumption).
     */
    const int ndf() const;
    const int constraints() const;
    /// normalized chi2 of the kinematic fit
    const float normalizedChi2() const;

    const ACFitParticle * tau() const;
    const ACFitParticle * prong0() const;
    const ACFitParticle * prong1() const;
    const ACFitParticle * prong2() const;
    const ACFitParticle * nu() const;

    /**
     returns summed lorentz vector of all visible daughters. this is the a1 for this specific decay assumption. per default the fitted state is regarded.
     */
    const TLorentzVector a1Lorentz(bool prefit = false) const;
    /// returns the GJ angle. per default the fitted state is regarded.
    const double gjAngle(bool prefit = false) const;
    /// returns the maximal allowed GJ angle. per default the fitted state is regarded.
    const double maxGJAngle(bool prefit = false) const;
    
    /// quality criterion: vertex significance between the rotated primary vertex and the secondary vertex of the tau decay
    const double vtxSignPVRotSV() const { return vtxSignPVRotSV_; }
    /// quality criterion: vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
    const double vtxSignPVRotPVRed() const { return vtxSignPVRotPVRed_; }
    /// quality criterion: mass of the a1 system
    const double a1Mass() const { return a1Mass_; }
    /// quality criterion: chi2 probability of the fit
    const double chi2Probability() const;
    /// size of tracks in the signal cone of the initial PFTau candidate
    const int sgnlConeTrkSize() const;
    /// quality criterion: transversal energy fraction between the intial PFTau and the final kinematic tau
    const double energyTFraction() const { return energyTFraction_; }


    ClassDef(ACFittedThreeProngDecay,1);

protected:
    /// constraints of the kinematic fit (called ndf in the KinematicFit package)
    int constraints_;
    /// real ndf depending on the decay mode. (WARNING: This is NOT the ndf() in RecoVertex/KinematicFit)
    int ndf_;
    
    /// quality criterion (may depend on decay mode): vertex significance between the rotated primary vertex and the secondary vertex of the tau decay (the tau carries rotated and reduced primVtx as initial vtx, but conversion into reco::Vertex needed)
    double vtxSignPVRotSV_;
    /// quality criterion (may depend on decay mode): vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
    double vtxSignPVRotPVRed_;
    /// quality criterion (may depend on decay mode): mass of the a1 system
    double a1Mass_;
    /// quality criterion (may depend on decay mode): transversal energy fraction between the intial PFTau and the final kinematic tau
    double energyTFraction_;

private:
    /// integrity check for correct amount of particles
    bool testDecay() const;
};

typedef ACFittedThreeProngDecay* pACFittedThreeProngDecay;

#endif
