#include "../interface/ACGenParticleConverter.h"

ACGenParticleConverter::ACGenParticleConverter(const reco::GenParticle& particle) {
    status_ = particle.status(); //is this meaningful?
    pdgId_ = particle.pdgId();
    charge_ = particle.charge();
    p4_ = TLorentzVector(particle.p4().Px(), particle.p4().Py(), particle.p4().Pz(), particle.p4().E());
    vertex_ = TVector3(particle.vertex().X(), particle.vertex().Y(), particle.vertex().Z());

    name_ = "other";
    if (abs(pdgId_)==25) name_ = "higgs";
    if (abs(pdgId_)==24) name_ = "W";
    if (abs(pdgId_)==23) name_ = "Z";
    if (abs(pdgId_)==15) name_ = "tau";
    if (abs(pdgId_)==16) name_ = "neutrino";
    if (abs(pdgId_)==11) {
        name_ = "electron";
        const reco::Candidate *genMother = particle.mother();
        if (genMother) {
            if (abs(genMother->pdgId())==22) name_ += "FromPhoton";
            else if (abs(genMother->pdgId())==111) name_ += "FromPi0";
            else if (abs(genMother->pdgId())!=15) {
//                printf("ACGenParticleConverter::ACGenParticleConverter: Rare electron mother %i\n", genMother->pdgId());
            }
        } else {
            printf("ACGenParticleConverter::ACGenParticleConverter: Bad link to mother!\n");
            throw 1;
        }
    }
    if (abs(pdgId_)==12) name_ = "electronnu";
    if (abs(pdgId_)==13) name_ = "muon";
    if (abs(pdgId_)==14) name_ = "muonnu";
    if (abs(pdgId_)==211) name_ = "pion";
    if (abs(pdgId_)==22) {
        name_ = "photon";
        const reco::Candidate *genMother = particle.mother();
        if (genMother) {
            if (abs(genMother->pdgId())==111) name_ += "FromPi0";
        } else {
            printf("ACGenParticleConverter::ACGenParticleConverter: Bad link to mother!\n");
            throw 1;
        }
    }
    if (abs(pdgId_) >= 1 && abs(pdgId_) <= 6) name_ = "quark";
    if (name_ == "other") {
        if (charge_==0) name_ = "neutral";
        else name_ = "charged";
    }
}

void ACGenParticleConverter::setGenDecayRef(const ACGenDecayRef & genDecayRef) {
    genDecayRef_ = genDecayRef;
}
