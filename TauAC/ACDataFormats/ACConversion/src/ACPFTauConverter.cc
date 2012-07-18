#include "../interface/ACPFTauConverter.h"

ACPFTauConverter::ACPFTauConverter(const edm::Event& evt, const reco::PFTauRef& tau, const std::vector<edm::InputTag> & discriminators, const PFJetConversionLog * conversionLogPFJet) {
    p4_ = TLorentzVector(tau->px(), tau->py(), tau->pz(), tau->energy());
    charge_ = tau->charge();
    for (std::vector<edm::InputTag>::const_iterator discr = discriminators.begin(); discr != discriminators.end(); ++discr) {
        edm::Handle<reco::PFTauDiscriminator> tmpHandle;
        evt.getByLabel(*discr, tmpHandle);
        if (!tmpHandle.isValid()) {
            edm::LogError("FinalTreeFiller") << "Invalid handle found in tau discriminators! Skip it.";
            continue;
        }
//        // make sure that PFTau and discriminators match each other. the tau ID might have been rerun, thus multiple copies of discriminators are in the event
//        const std::string & process = discr->provenance()->processName();
//        const edm::Provenance & prov = evt.getProvenance(tau.id());
//        if (process != prov.processName()) continue;
//        const std::string & label = discr->provenance()->moduleLabel();
//        if (label.substr(0, 24) != pfTauDiscriminatorPattern) continue;
//        //std::cout<<"Discriminator "<<label<<std::endl;
        try {
            discriminators_.insert(std::make_pair(discr->label(), (*tmpHandle)[tau]));
            //discriminators_.insert(std::make_pair(label, (**discr)[tau]));
        } catch (edm::errors::ErrorCodes error) {
            std::cout<<"ACPFTauConverter::ACPFTauConverter: Catched an error "<<error<<" (expecting "<<edm::errors::InvalidReference<<")"<<std::endl;
        } catch (...) {
            std::cout<<"ACPFTauConverter::ACPFTauConverter: Catched an unknown error (expecting "<<edm::errors::InvalidReference<<")"<<std::endl;            
        }
    }
    //std::cout<<"stored "<<discriminators_.size()<<" discriminators"<<std::endl;

    const ACJetRef * jetRef = 0;
    if (tau->jetRef().isAvailable()) {
        jetRef = conversionLogPFJet->getConverted(tau->jetRef());
    }
    // reference might be missing due to candidate filters
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
