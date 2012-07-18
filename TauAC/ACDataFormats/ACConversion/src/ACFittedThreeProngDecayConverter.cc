#include "../interface/ACFittedThreeProngDecayConverter.h"

ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter(const edm::Event& evt, const SelectedKinematicDecay & decay, MCTruthMatching * kinematicParticleMatching, std::vector<ACGenDecay *> * genTauDecays, std::vector<ACFitParticle *> * kinematicParticles, PFTauConversionLog * conversionLogPFTau) {
    //kinematicParticleMatching is deprecated!!! remove soon, only its counters are still used
    
    iterations_ = decay.iterations();
    csum_ = decay.csum();
    maxiterations_ = decay.maxiterations();
    mincsum_ = decay.mincsum();
    chi2_ = decay.chi2();
    constraints_ = decay.constraints();
    ndf_ = decay.ndf();

    const ACPFTauRef * pfTauRef = conversionLogPFTau->getConverted(decay.PFTauRef());
    // reference might be missing
    if (pfTauRef == 0) {
        PFTauRef_ = ACPFTauRef();
    } else {
        PFTauRef_ = *pfTauRef;
    }

    SelectedKinematicParticleCollection decayparticles = decay.particles();
    if (decayparticles.size() > 0) {
		bool motherMissmatch = false;
		const ACGenDecayRef * matchedGenDecay = 0;
        for (SelectedKinematicParticleCollection::const_iterator particle=decayparticles.begin(); particle!=decayparticles.end(); ++particle) {
            ACGenParticleRef genRef;
            /// try to match all charged daughters of a tau decay
            if (particle->charge() != 0 && particle != decayparticles.begin()) {
                //std::cout<<"Work on pion candidate with charge "<<particle->charge()<<" and pt "<<particle->candRef()->pt()<<":"<<std::endl;
                //genRef = kinematicParticleMatching->getMatching(evt, &*particle);//deprecated as only the closest track is considered (poor efficiency on pileup)
                // instead: assign all objects by closest match in dR that also matches in charge and (roughly) in pt
                // WARNING: ambigious assignments possible
                for (std::vector<ACGenDecay *>::const_iterator genTauDecay = genTauDecays->begin(); genTauDecay != genTauDecays->end(); ++genTauDecay) {
                    std::vector<const ACGenParticle*> chargedGenDaughters;
                    (*genTauDecay)->chargedDaughters(chargedGenDaughters);
                    for (std::vector<const ACGenParticle *>::const_iterator genParticle = chargedGenDaughters.begin(); genParticle != chargedGenDaughters.end(); ++genParticle) {
                        if (particle->charge()!=(*genParticle)->charge()) continue;
//                        if (fabs(particle->candRef()->pt() - (*genParticle)->pt()) > .2 * (*genParticle)->pt()) continue;
                        const double & dR = reco::deltaR(*(particle->candRef()), **genParticle);
                        if (dR > 0.01) continue;
                        //std::cout<<"\tFound a similar track from pdgID "<<(*genParticle)->pdgId()<<" at dR = "<<dR<<" with pT = "<<(*genParticle)->pt()<<std::endl;
                        if (genRef.isValid()) {
                            if (genRef.get() != *genParticle) {
                                const double & dRInitial = reco::deltaR(*(particle->candRef()), *genRef);
                                if (dR < dRInitial) {
                                    //std::cout<<"\t\tReplace initial match to "<<genRef->pdgId()<<" at dR = "<<dRInitial<<std::endl;
                                    //update matching now
                                    genRef = ACGenParticleRef(*genParticle);
                                } else {
                                    //std::cout<<"\t\tKeep initial match to "<<genRef->pdgId()<<" at dR = "<<dRInitial<<std::endl;
                                }
                            }
                        } else {
                            //std::cout<<"\t\tStore match."<<std::endl;
                            //update matching now
                            genRef = ACGenParticleRef(*genParticle);
                        }
                    }
                }
                kinematicParticleMatching->increment(genRef.isValid());//just a counter
                //if (!genRef.isValid()) std::cout<<"\tNo match found!"<<std::endl;

                //now scan for common mother
                if (genRef.isValid()) {
                    if (!motherMissmatch) {
                        if (!matchedGenDecay) {
                            matchedGenDecay = &(genRef->genDecayRef());
                            if (!matchedGenDecay->isValid()) {
                                printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: Error! Bad generator-decay reference!\n");
                                motherMissmatch = true;
                            }
                        } else {
                            if (matchedGenDecay->get() != genRef->genDecayRef().get()) {
                                printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: WARNING! decays missmatched. gentau %s != gentau %s\n", genRef->genDecayRef()->name().c_str(), (*matchedGenDecay)->name().c_str());
                                matchedGenDecay = 0;
                                motherMissmatch = true;
                            } else {
                                //printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: decays matched. gentau %f == gentau %f\n", genRef->genDecayRef()->topParticle()->pt(), (*matchedGenDecay)->topParticle()->pt());
                            }
                        }
                    }
                }
            }
            ACFitParticleConverter tmp(*particle, decayparticles.begin()->charge(), &genRef);
            ACFitParticle * tmpP = new ACFitParticle();
            *tmpP = tmp;
            /// store particles of a kinematic fit into a branch
            kinematicParticles->push_back(tmpP);
            /// store reference to particle into decay
            particles_.push_back(ACFitParticleRef(kinematicParticles->back()));
        }
        // Accept unambigious assignment also if decay modes mismatch. I.e., consider a decay as matched if at least one charged particle was matched.
        if (matchedGenDecay != 0) {
            // also store decays with a wrong charge
            genDecayRef_ = *matchedGenDecay;
            // test for same charge of the top particle
//            if (particles()->front()->charge() == (*matchedGenDecay)->particles()->front()->charge()) {
//                printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: common tau found with charge %f, pt %f\n", (*matchedGenDecay)->particles()->front()->charge(), (*matchedGenDecay)->particles()->front()->pt());
//            } else printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: WARNING! Store matching with bad tau charge!\n");
		} //else printf("ACFittedThreeProngDecayConverter::ACFittedThreeProngDecayConverter: no common decay found!\n");
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
    
    /// store quality discriminators that cannot directly be calculated from stored members only (e.g. conversion into reco::Vertex format would be needed). FIXME: replace this by  a dynamic calculation (depending on the decay mode)
    vtxSignPVRotSV_ = decay.vtxSignPVRotSV();
    vtxSignPVRotPVRed_ = decay.vtxSignPVRotPVRed();
    a1Mass_ = decay.a1Mass();
    energyTFraction_ = decay.energyTFraction();
}
