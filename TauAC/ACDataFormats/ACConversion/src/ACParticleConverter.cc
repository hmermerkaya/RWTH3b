#include "../interface/ACParticleConverter.h"

ACParticleConverter::ACParticleConverter(const reco::Muon& particle, ACGenParticleRef genRef) {
    status_ = particle.status(); //is this meaningful?
    name_ = "muon";
    if (particle.charge() < 0) {
        pdgId_ = -13;
    } else {
        pdgId_ = 13;
    }
    charge_ = particle.charge();
    p4_ = TLorentzVector(particle.p4().Px(), particle.p4().Py(), particle.p4().Pz(), particle.p4().E());
    vertex_ = TVector3(particle.vertex().X(), particle.vertex().Y(), particle.vertex().Z());
    matrix_.ResizeTo(TMatrixDSym(7)); //empty matrix at the moment
    genRef_ = genRef;
    /// to be ignored:    ACGenDecayRef genDecayRef_
}
ACParticleConverter::ACParticleConverter(const reco::GsfElectron& particle, ACGenParticleRef genRef) {
    status_ = particle.status(); //is this meaningful?
    name_ = "electron";
    if (particle.charge() < 0) {
        pdgId_ = -11;
    } else {
        pdgId_ = 11;
    }
    charge_ = particle.charge();
    p4_ = TLorentzVector(particle.p4().Px(), particle.p4().Py(), particle.p4().Pz(), particle.p4().E());
    vertex_ = TVector3(particle.vertex().X(), particle.vertex().Y(), particle.vertex().Z());
    matrix_.ResizeTo(TMatrixDSym(7)); //empty matrix at the moment
    genRef_ = genRef;
    /// to be ignored: ACGenDecayRef genDecayRef_
}