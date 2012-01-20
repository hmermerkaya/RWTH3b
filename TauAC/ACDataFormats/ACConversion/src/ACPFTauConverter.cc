#include "../interface/ACPFTauConverter.h"

ACPFTauConverter::ACPFTauConverter(const edm::Event& evt, const reco::PFTauRef& tau, const std::vector<edm::InputTag> & tauDiscriminators, const PFJetConversionLog * conversionLogPFJet) {
    p4_ = TLorentzVector(tau->px(), tau->py(), tau->pz(), tau->energy());
    charge_ = tau->charge();
    for (std::vector<edm::InputTag>::const_iterator discr = tauDiscriminators.begin(); discr != tauDiscriminators.end(); ++discr) {
        edm::Handle<reco::PFTauDiscriminator> tmpHandle;
        evt.getByLabel(*discr, tmpHandle);
        
        discriminators_.insert(std::make_pair(discr->label(), (*tmpHandle)[tau]));
    }

    const ACJetRef * jetRef = 0;
    if (tau->jetRef().isAvailable()) {
        jetRef = conversionLogPFJet->getConverted(tau->jetRef());
    }
    // reference might be missing
    if (jetRef == 0) {
        jetRef_ = ACJetRef();
    } else {
        jetRef_ = *jetRef;
    }

    signalPFChargedHadrCands_ = tau->signalPFChargedHadrCands().size();
    signalPFNeutrHadrCands_ = tau->signalPFNeutrHadrCands().size();
    signalPFGammaCands_ = tau->signalPFGammaCands().size();
    signalPiZeroCands_ = tau->signalPiZeroCandidates().size();

    isolationPFChargedHadrCands_ = tau->isolationPFChargedHadrCands().size();
    isolationPFNeutrHadrCands_ = tau->isolationPFNeutrHadrCands().size();
    isolationPFGammaCands_ = tau->isolationPFGammaCands().size();
    isolationPiZeroCands_ = tau->isolationPiZeroCandidates().size();
    
}
