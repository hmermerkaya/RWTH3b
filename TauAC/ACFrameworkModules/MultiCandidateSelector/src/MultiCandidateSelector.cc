#include "ACFrameworkModules/MultiCandidateSelector/interface/MultiCandidateSelector.h"

MultiCandidateSelector::MultiCandidateSelector(const edm::ParameterSet& iConfig) {
    filters_ = iConfig.getParameter<std::vector<edm::ParameterSet> >("filters");
    
    registerProducts();
}
MultiCandidateSelector::~MultiCandidateSelector() {}

void MultiCandidateSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    for (std::vector<edm::ParameterSet>::const_iterator iter = filters_.begin(); iter != filters_.end(); ++iter) {
        std::string thisType = iter->getParameter<std::string>("type");
        edm::InputTag thisInputTag = iter->getParameter<edm::InputTag>("src");
        std::string thisCut = iter->getParameter<std::string>("cut");
        
        if (thisType == "muons") {
            edm::Handle<reco::MuonCollection> pIn;
            iEvent.getByLabel(thisInputTag, pIn);
            
            std::auto_ptr<reco::MuonCollection> pOut(new reco::MuonCollection);
            StringCutObjectSelector<reco::Muon> thisSelector(thisCut);
            
            for (reco::MuonCollection::const_iterator i = pIn->begin(); i != pIn->end(); ++i) {
                if (thisSelector(*i)) pOut->push_back(*i);
            }
            iEvent.put(pOut, thisInputTag.label());
        } else if (thisType == "electrons") {
            edm::Handle<reco::GsfElectronCollection> pIn;
            iEvent.getByLabel(thisInputTag, pIn);
            
            std::auto_ptr<reco::GsfElectronCollection> pOut(new reco::GsfElectronCollection);
            StringCutObjectSelector<reco::GsfElectron> thisSelector(thisCut);
            
            for (reco::GsfElectronCollection::const_iterator i = pIn->begin(); i != pIn->end(); ++i) {
                if (thisSelector(*i)) pOut->push_back(*i);
            }
            iEvent.put(pOut, thisInputTag.label());
        } else if (thisType == "pfTaus") {
            edm::Handle<reco::PFTauCollection> pIn;
            iEvent.getByLabel(thisInputTag, pIn);
            
            std::auto_ptr<reco::PFTauCollection> pOut(new reco::PFTauCollection);
            StringCutObjectSelector<reco::PFTau> thisSelector(thisCut);
            
            for (reco::PFTauCollection::const_iterator i = pIn->begin(); i != pIn->end(); ++i) {
                if (thisSelector(*i)) pOut->push_back(*i);
            }
            iEvent.put(pOut, thisInputTag.label());
        } else if (thisType == "pfJets") {
            edm::Handle<reco::PFJetCollection> pIn;
            iEvent.getByLabel(thisInputTag, pIn);
            
            std::auto_ptr<reco::PFJetCollection> pOut(new reco::PFJetCollection);
            StringCutObjectSelector<reco::PFJet> thisSelector(thisCut);
            
            for (reco::PFJetCollection::const_iterator i = pIn->begin(); i != pIn->end(); ++i) {
                if (thisSelector(*i)) pOut->push_back(*i);
            }
            iEvent.put(pOut, thisInputTag.label());
        } else if (thisType == "tracks") {
            edm::LogWarning("MultiCandidateSelector") << "Type '" << thisType << "' not implemented yet!";
        } else {
            edm::LogWarning("MultiCandidateSelector") << "Type '" << thisType << "' not known!";
        }
    }
}
void MultiCandidateSelector::registerProducts() {
    for (std::vector<edm::ParameterSet>::const_iterator iter = filters_.begin(); iter != filters_.end(); ++iter) {
        std::string thisType = iter->getParameter<std::string>("type");
        edm::InputTag thisInputTag = iter->getParameter<edm::InputTag>("src");
        
        if (thisType == "muons") {
            produces<reco::MuonCollection>(thisInputTag.label());
        } else if (thisType == "electrons") {
            produces<reco::GsfElectronCollection>(thisInputTag.label());
        } else if (thisType == "pfTaus") {
            produces<reco::PFTauCollection>(thisInputTag.label());
        } else if (thisType == "pfJets") {
            produces<reco::PFJetCollection>(thisInputTag.label());
        } else if (thisType == "tracks") {
            produces<reco::TrackCollection>(thisInputTag.label());
        } else {
            edm::LogWarning("MultiCandidateSelector") << "Type '" << thisType << "' not known!";
        }
    }
}
void MultiCandidateSelector::beginJob() {
}
void MultiCandidateSelector::endJob() {
}
void MultiCandidateSelector::beginRun(edm::Run&, edm::EventSetup const&) {
}
void MultiCandidateSelector::endRun(edm::Run&, edm::EventSetup const&) {
}
void MultiCandidateSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
void MultiCandidateSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
void MultiCandidateSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    //Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MultiCandidateSelector);
