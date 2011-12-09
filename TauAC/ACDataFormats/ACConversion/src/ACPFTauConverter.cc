#include "../interface/ACPFTauConverter.h"

ACPFTauConverter::ACPFTauConverter(const edm::Event& evt, const reco::PFTauRef& tau, const std::vector<edm::InputTag> & tauDiscriminators) {
    p4_ = TLorentzVector(tau->px(), tau->py(), tau->pz(), tau->energy());
    charge_ = tau->charge();
    signalPFChargedHadrCands_ = tau->signalPFChargedHadrCands().size();
    signalPFNeutrHadrCands_ = tau->signalPFNeutrHadrCands().size();

    for (std::vector<edm::InputTag>::const_iterator discr = tauDiscriminators.begin(); discr != tauDiscriminators.end(); ++discr) {
        edm::Handle<reco::PFTauDiscriminator> tmpHandle;
        evt.getByLabel(*discr, tmpHandle);

        discriminators_.insert(std::make_pair(discr->label(), (*tmpHandle)[tau]));
    }
}
