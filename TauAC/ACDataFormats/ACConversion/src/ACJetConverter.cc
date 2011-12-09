#include "../interface/ACJetConverter.h"

ACJetConverter::ACJetConverter(const reco::Jet& jet) {
    p4_ = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy());
    charge_ = jet.charge();
    multiplicity_ = jet.nConstituents();
}
ACJetConverter::ACJetConverter(const reco::PFJet& jet) {
    p4_ = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy());
    charge_ = jet.charge();
    multiplicity_ = jet.nConstituents();
}