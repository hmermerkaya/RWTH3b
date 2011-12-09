#include "../interface/ACEvent.h"

ACEvent::ACEvent() {
    /// pointers to branch content have to be set to 0 to prevent from root crash
    eventInfo_ = 0;
    eventGlobals_ = 0;
    trigger_ = 0;
    offlinePV_ = 0;
    reducedPV_ = 0;
    generator_ = 0;
    genTauDecays_ = 0;
    muons_ = 0;
    electrons_ = 0;
    tauDecays_ = 0;
    fittedThreeProngParticles_ = 0;
    pfJets_ = 0;
    pfTaus_ = 0;
}

const ACEventInfo * ACEvent::eventInfo() const {
    checkContent(eventInfo_, "eventInfo");
    return eventInfo_;
}
ACEventInfo ** ACEvent::linkEventInfo() {
    return &eventInfo_;
}

const ACEventGlobals * ACEvent::eventGlobals() const {
    checkContent(eventGlobals_, "eventGlobals");
    return eventGlobals_;
}
ACEventGlobals ** ACEvent::linkEventGlobals() {
    return &eventGlobals_;
}

const ACTrigger * ACEvent::trigger() const {
    checkContent(trigger_, "trigger");
    return trigger_;
}
ACTrigger ** ACEvent::linkTrigger() {
    return &trigger_;
}

const std::vector<ACVertex *>& ACEvent::offlinePV() const {
    checkContent(offlinePV_, "offlinePV");
    return *offlinePV_;
}
std::vector<ACVertex *> ** ACEvent::linkOfflinePV() {
    return &offlinePV_;
}

const std::vector<ACVertex *>& ACEvent::reducedPV() const {
    checkContent(reducedPV_, "reducedPV");
    return *reducedPV_;
}
std::vector<ACVertex *> ** ACEvent::linkReducedPV() {
    return &reducedPV_;
}

const std::vector<ACGenParticle *>& ACEvent::generator() const {
    checkContent(generator_, "generator");
    return *generator_;
}
std::vector<ACGenParticle *> ** ACEvent::linkGenerator() {
    return &generator_;
}

const std::vector<ACGenDecay *>& ACEvent::genTauDecays() const {
    checkContent(genTauDecays_, "genTauDecays");
    return *genTauDecays_;
}
std::vector<ACGenDecay *> ** ACEvent::linkGenTauDecays() {
    return &genTauDecays_;
}

const std::vector<ACParticle *>& ACEvent::muons() const {
    checkContent(muons_, "muons");
    return *muons_;
}
std::vector<ACParticle *> ** ACEvent::linkMuons() {
    return &muons_;
}

const std::vector<ACParticle *>& ACEvent::electrons() const {
    checkContent(electrons_, "electrons");
    return *electrons_;
}
std::vector<ACParticle *> ** ACEvent::linkElectrons() {
    return &electrons_;
}

const std::vector<ACFittedThreeProngDecay *>& ACEvent::tauDecays() const {
    checkContent(tauDecays_, "tauDecays");
    return *tauDecays_;
}
std::vector<ACFittedThreeProngDecay *> ** ACEvent::linkTauDecays() {
    return &tauDecays_;
}

std::vector<ACFitParticle *> ** ACEvent::linkFittedThreeProngParticles() {
    return &fittedThreeProngParticles_;
}

const std::vector<ACJet *>& ACEvent::pfJets() const {
    checkContent(pfJets_, "pfJets");
    return *pfJets_;
}
std::vector<ACJet *> ** ACEvent::linkPFJets() {
    return &pfJets_;
}

const std::vector<ACPFTau *>& ACEvent::pfTaus() const {
    checkContent(pfTaus_, "pfTaus");
    return *pfTaus_;
}
std::vector<ACPFTau *> ** ACEvent::linkPFTaus() {
    return &pfTaus_;
}
