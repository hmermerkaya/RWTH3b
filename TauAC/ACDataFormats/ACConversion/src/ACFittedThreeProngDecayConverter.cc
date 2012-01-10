#include "../interface/ACFittedThreeProngDecayConverter.h"

ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter(const edm::Event& evt, const SelectedKinematicDecay & decay, MCTruthMatching * kinematicParticleMatching, std::vector<ACFitParticle *> * kinematicParticles, PFTauMatching * pfTauMatching) {
    if (!decay.topParticle()) {
        printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: ERROR! Invalid top particle! Cannot call the constructor.\n");
        return;
    }
    iterations_ = decay.topParticle()->iterations();
    csum_ = decay.topParticle()->csum();
    maxiterations_ = decay.topParticle()->maxiterations();
    mincsum_ = decay.topParticle()->mincsum();
    chi2_ = decay.topParticle()->chi2();
    constraints_ = decay.topParticle()->ndf();

    const ACPFTauRef * pfTauRef = 0;
    if (decay.PFTauRef().size() > 0) {
        pfTauRef = pfTauMatching->getConverted(decay.PFTauRef().front());
    }
    if (pfTauRef == 0) {
        PFTauRef_ = ACPFTauRef();
    } else {
        PFTauRef_ = *pfTauRef;
    }

    std::vector<const SelectedKinematicParticle *> decayparticles;
    decay.particles(decayparticles);
    if (decayparticles.size() > 0) {
		bool motherMissmatch = false;
		const ACGenDecayRef * matchedGenDecay = 0;
        for (std::vector<const SelectedKinematicParticle *>::const_iterator particle=decayparticles.begin(); particle!=decayparticles.end(); ++particle) {
            const ACGenParticleRef * genRef = 0;
            /// try to match all charged daughters of a tau decay
            if ((*particle)->charge() != 0 && particle != decayparticles.begin()) {
                genRef = kinematicParticleMatching->getMatching(evt, *particle);
                if (!genRef) matchedGenDecay = 0;
                else {
                    if (!motherMissmatch) {
                        if (matchedGenDecay==0) {
                            matchedGenDecay = &((*genRef)->genDecayRef());
                            if (!matchedGenDecay->isValid()) {
                                printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: Error! Bad generator-decay reference!\n");
                                motherMissmatch = true;
                            }
                        } else {
                            if (matchedGenDecay != &((*genRef)->genDecayRef())) {
                                //printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: decays missmatched. gentau %s != gentau %s\n", (*genRef)->genDecayRef()->name().c_str(), (*matchedGenDecay)->name().c_str());
                                matchedGenDecay = 0;
                                motherMissmatch = true;
                            }
                        }
                    }
                }
            }
            ACFitParticleConverter tmp(**particle, (*(decayparticles.begin()))->charge(), genRef);
            ACFitParticle * tmpP = new ACFitParticle();
            *tmpP = tmp;
            /// store particles of a kinematic fit into a branch
            kinematicParticles->push_back(tmpP);
            /// store reference to particle into decay
            particles_.push_back(ACFitParticleRef(kinematicParticles->back()));
        }
        // Accept unambigious assignment also if decay modes mismatch. I.e., consider a decay as matched if at least one charged particle was matched.
        if (matchedGenDecay != NULL) {
            // test for same charge of the top particle
            if (particles()->front()->charge() == (*matchedGenDecay)->particles()->front()->charge()) {
                genDecayRef_ = *matchedGenDecay;
            } //else printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: bad charge!\n");
		}// else printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: no common decay found!\n");
    } else {
        printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: ERROR! No particles in SelectedKinematicDecay.\n");
    }
    if (particles()->size() > 0) {
        if (particles()->front()->charge() == -1) {
            name_ = std::string("Tau-");
        } else {
            name_ = std::string("Tau+");
        }
    }
}
