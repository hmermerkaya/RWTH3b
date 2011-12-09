#include "ACFrameworkModules/FinalTreeFiller/interface/FinalTreeFiller.h"

FinalTreeFiller::FinalTreeFiller(const edm::ParameterSet& iConfig):
genSignalTag_(iConfig.getParameter<edm::InputTag>("genSignal")),
chargedTauDaughterMatchMapTag_(iConfig.getParameter<edm::InputTag>("chargedTauDaughterMatchMap")),
primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
reducedPrimVtxTag_(iConfig.getParameter<edm::InputTag>("reducedPrimVtx")),
pileupInfoTag_(iConfig.getParameter<edm::InputTag>("pileupInfo")),
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
muonTag_(iConfig.getParameter<edm::InputTag>("muons")),
electronTag_(iConfig.getParameter<edm::InputTag>("electrons")),
kinematicTausTag_(iConfig.getParameter<edm::InputTag>("kinematicTaus")),
pfMETTag_(iConfig.getParameter<edm::InputTag>("pfMET")),
tcMETTag_(iConfig.getParameter<edm::InputTag>("tcMET")),
pfJetTag_(iConfig.getParameter<edm::InputTag>("pfJets")),
pfTauTag_(iConfig.getParameter<edm::InputTag>("pfTaus")),
pfTauDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("pfTauDiscriminators"))
{
    hltChanged_ = true;
    eventTree_ = 0;
    acEventInfo_ = 0;
    acEventGlobals_ = 0;
    acTrigger_ = 0;
    offlinePV_ = 0;
    reducedPV_ = 0;
    muons_ = 0;
    electrons_ = 0;
    kinematicParticles_ = 0;
    kinematicDecays_ = 0;
    genParticles_ = 0;
    genTauDecays_ = 0;
    pfJets_ = 0;
    pfTaus_ = 0;
}
FinalTreeFiller::~FinalTreeFiller() {}

void FinalTreeFiller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    evtCnt_++;

    /// get current TObject count
    unsigned int objectCount = TProcessID::GetObjectCount();

    storeEvent(iEvent);
    eventTree_->Fill();

    /// reset TObject count to avoid slow downs
    TProcessID::SetObjectCount(objectCount);

    /// clear memory
    deleteVectorOfPointers(offlinePV_);
    deleteVectorOfPointers(reducedPV_);
    deleteVectorOfPointers(muons_);
    deleteVectorOfPointers(electrons_);
    deleteVectorOfPointers(kinematicParticles_);
    deleteVectorOfPointers(kinematicDecays_);
    deleteVectorOfPointers(genParticles_);
    deleteVectorOfPointers(genTauDecays_);
    deleteVectorOfPointers(pfJets_);
    deleteVectorOfPointers(pfTaus_);
}
void FinalTreeFiller::beginJob() {
    evtCnt_ = 0;
    runCnt_ = 0;

    kinematicParticleMatching_ = new MCTruthMatching(chargedTauDaughterMatchMapTag_);
    pfTauMatching_ = new PFTauMatching();

    eventTree_ = fileService_->make<TTree>("TauACEvent", "TauAC tree filled each event");
    //runTree_ = fileService_->make<TTree>("TauACRun", "TauAC tree filled each run");

    acEventInfo_ = new ACEventInfo();
    eventTree_->Bronch("ACEventInfo", "ACEventInfo", &acEventInfo_, 32000, 0);

    acEventGlobals_ = new ACEventGlobals();
    eventTree_->Bronch("ACEventGlobals", "ACEventGlobals", &acEventGlobals_, 32000, 0);

    acTrigger_ = new ACTrigger();
    eventTree_->Bronch("ACTrigger", "ACTrigger", &acTrigger_, 32000, 0);

    offlinePV_ = new std::vector<ACVertex *>();
    eventTree_->Bronch("offlinePV", "std::vector<ACVertex *>", &offlinePV_, 32000, 0);
    reducedPV_ = new std::vector<ACVertex *>();
    eventTree_->Bronch("reducedPV", "std::vector<ACVertex *>", &reducedPV_, 32000, 0);

    genParticles_ = new std::vector<ACGenParticle *>();
    eventTree_->Bronch("generator", "std::vector<ACGenParticle *>", &genParticles_, 32000, 0);

    genTauDecays_ = new std::vector<ACGenDecay *>();
    eventTree_->Bronch("genTauDecays", "std::vector<ACGenDecay *>", &genTauDecays_, 32000, 0);

    muons_ = new std::vector<ACParticle *>();
    eventTree_->Bronch("Muons", "std::vector<ACParticle *>", &muons_, 32000, 0);

    electrons_ = new std::vector<ACParticle *>();
    eventTree_->Bronch("Electrons", "std::vector<ACParticle *>", &electrons_, 32000, 0);

    kinematicParticles_ = new std::vector<ACFitParticle *>();
    eventTree_->Bronch("ACFittedThreeProngParticles", "std::vector<ACFitParticle *>", &kinematicParticles_, 32000, 0);

    kinematicDecays_ = new std::vector<ACFittedThreeProngDecay *>();
    eventTree_->Bronch("ACFittedThreeProngDecays", "std::vector<ACFittedThreeProngDecay *>", &kinematicDecays_, 32000, 0);

    pfJets_ = new std::vector<ACJet *>();
    eventTree_->Bronch("PFJets", "std::vector<ACJet *>", &pfJets_, 32000, 0);

    pfTaus_ = new std::vector<ACPFTau *>();
    eventTree_->Bronch("PFTaus", "std::vector<ACPFTau *>", &pfTaus_, 32000, 0);
}
void FinalTreeFiller::endJob() {
    kinematicParticleMatching_->printOutro();

    /// clear memory
    delete kinematicParticleMatching_;
    delete acEventInfo_;
    delete acEventGlobals_;
    delete acTrigger_;
    delete offlinePV_;
    delete reducedPV_;
    delete muons_;
    delete electrons_;
    delete kinematicParticles_;
    delete kinematicDecays_;
    delete genParticles_;
    delete genTauDecays_;
    delete pfJets_;
    delete pfTaus_;
}
void FinalTreeFiller::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    /// initialize HLTConfigProvider
    HLTCP_.init(iRun,iSetup,triggerResultsTag_.process(),hltChanged_);
    runCnt_++;
}
void FinalTreeFiller::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

}
void FinalTreeFiller::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {

}
void FinalTreeFiller::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {

}
void FinalTreeFiller::storeEvent(const edm::Event& evt) {
    kinematicParticleMatching_->nextEvt();
    pfTauMatching_->nextEvt();

    /// set event information
    edm::Handle<reco::GenParticleCollection> genStatus;
    /// check whether generator information is present in the event
    if (evt.getByLabel(edm::InputTag("genParticles"), genStatus)) {
        *acEventInfo_ = ACEventInfo(evt.id().event(), evt.id().luminosityBlock(), evt.id().run(), "MC");
    } else {
        *acEventInfo_ = ACEventInfo(evt.id().event(), evt.id().luminosityBlock(), evt.id().run(), "data");
    }

    /// set event globals
    edm::Handle<reco::PFMETCollection> pfMET;
    edm::Handle<reco::METCollection> tcMET;
    if (loadCollection(evt, pfMETTag_, pfMET) && loadCollection(evt, tcMETTag_, tcMET)) {
        *acEventGlobals_ = ACEventGlobals(TVector3(pfMET->begin()->px(), pfMET->begin()->py(), 0.0), TVector3(tcMET->begin()->px(), tcMET->begin()->py(), 0.0), pfMET->begin()->sumEt(), tcMET->begin()->sumEt());
    }

    /// set trigger results
    edm::Handle<edm::TriggerResults> HLTR;
    if (loadCollection(evt, triggerResultsTag_, HLTR)) {
        if (HLTCP_.inited()) {
            /// generate trigger path map
            edm::TriggerNames triggerNames = evt.triggerNames(*HLTR);
            std::map<std::string,int> hltPathMap;
            for (unsigned int i = 0; i != triggerNames.triggerNames().size(); ++i) {
                hltPathMap[triggerNames.triggerName(i)] = HLTR->accept(i);
            }
            *acTrigger_ = ACTrigger(HLTCP_.tableName(), hltPathMap);
        }
    }

    /// set primary-vertex collection
    edm::Handle<reco::VertexCollection> offlinePrimaryVertices;
    if (loadCollection(evt, primVtxTag_, offlinePrimaryVertices)) {
        *offlinePV_ = std::vector<ACVertex *>();
        for (reco::VertexCollection::const_iterator candidate = offlinePrimaryVertices->begin(); candidate != offlinePrimaryVertices->end(); ++candidate) {
            ACVertexConverter tmp(*candidate);
            ACVertex * tmpP = new ACVertex();
            *tmpP = tmp;
            offlinePV_->push_back(tmpP);
        }
    }

    /// set primary-vertex collection (recalculated w/o the tracks belonging to kinematic taus)
    edm::Handle<reco::VertexCollection> reducedPrimaryVertices;
    if (loadCollection(evt, reducedPrimVtxTag_, reducedPrimaryVertices)) {
        *reducedPV_ = std::vector<ACVertex *>();
        for (reco::VertexCollection::const_iterator candidate = reducedPrimaryVertices->begin(); candidate != reducedPrimaryVertices->end(); ++candidate) {
            ACVertexConverter tmp(*candidate);
            ACVertex * tmpP = new ACVertex();
            *tmpP = tmp;
            reducedPV_->push_back(tmpP);
        }
    }

    /// set gerenator information
    edm::Handle<reco::GenParticleCollection> genCands;
    if (loadCollection(evt, genSignalTag_, genCands)) {
        *genParticles_ = std::vector<ACGenParticle *>();
        *genTauDecays_ = std::vector<ACGenDecay *>();
        ACGenDecayRef latestGenDecayRef;
        unsigned int index = 0;
        for (reco::GenParticleCollection::const_iterator candidate = genCands->begin(); candidate != genCands->end(); ++candidate, index++) {
            ACGenParticleConverter tmp(*candidate);
            ACGenParticle * tmpP = new ACGenParticle();
            *tmpP = tmp;
            genParticles_->push_back(tmpP);
            ACGenParticleRef tmpRef(genParticles_->back());
            kinematicParticleMatching_->logConversion(reco::GenParticleRef(genCands, index), tmpRef);
            /// store every occuring tau decay, expect the tau to be followed by its daughters
            if (abs(candidate->pdgId())==15) {
                ACGenDecayConverter tmpDecay(tmpRef);
                ACGenDecay * tmpDecayP = new ACGenDecay();
                *tmpDecayP = tmpDecay;
                genTauDecays_->push_back(tmpDecayP);
            } else {
                /// skip all particles before the first tau.
                if (genTauDecays_->size()<1) continue;
                /// skip quarks
                if (abs(candidate->pdgId()) >= 1 && abs(candidate->pdgId()) <= 6) continue;
                ACGenDecayConverter * tmpConverter = static_cast<ACGenDecayConverter*>(genTauDecays_->back());
                if (tmpConverter) {
                    tmpConverter->addDaughter(tmpRef);
                }
            }
            /// store the reference to the generator tau decay that contains the current particle
            ACGenParticleConverter * tmpConverter = static_cast<ACGenParticleConverter*>(genParticles_->back());
            if (tmpConverter) {
                tmpConverter->setGenDecayRef(ACGenDecayRef(genTauDecays_->back()));
            }
        }
    }

    /// set muon collection
    edm::Handle<reco::MuonCollection> muons;
    if (loadCollection(evt, muonTag_, muons)) {
        *muons_ = std::vector<ACParticle *>();
        for (reco::MuonCollection::const_iterator imu = muons->begin(); imu != muons->end(); ++imu) {
            ACParticleConverter tmp(*imu);
            ACParticle * tmpP = new ACParticle();
            *tmpP = tmp;
            muons_->push_back(tmpP);
        }
    }

    /// set electron collection
    edm::Handle<reco::GsfElectronCollection> electrons;
    if (loadCollection(evt, electronTag_, electrons)) {
        *electrons_ = std::vector<ACParticle *>();
        for (reco::GsfElectronCollection::const_iterator ie = electrons->begin(); ie != electrons->end(); ++ie) {
            ACParticleConverter tmp(*ie);
            ACParticle * tmpP = new ACParticle();
            *tmpP = tmp;
            electrons_->push_back(tmpP);
        }
    }

    /// set pf jet collection
    edm::Handle<reco::PFJetCollection> pfJets;
    if (loadCollection(evt, pfJetTag_, pfJets)) {
        *pfJets_ = std::vector<ACJet *>();
        for (reco::PFJetCollection::const_iterator ijet = pfJets->begin(); ijet != pfJets->end(); ++ijet) {
            ACJetConverter tmp(*ijet);
            ACJet * tmpJet = new ACJet();
            *tmpJet = tmp;
            pfJets_->push_back(tmpJet);
        }
    }

    /// set pf tau collection
    edm::Handle<reco::PFTauCollection> pfTaus;
    if (loadCollection(evt, pfTauTag_, pfTaus)) {
        *pfTaus_ = std::vector<ACPFTau *>();
        unsigned int index = 0;
        for (reco::PFTauCollection::const_iterator itau = pfTaus->begin(); itau != pfTaus->end(); ++itau, index++) {
            ACPFTauConverter tmp(evt, reco::PFTauRef(pfTaus, index), pfTauDiscriminatorTags_);
            ACPFTau * tmpTau = new ACPFTau();
            *tmpTau = tmp;
            pfTaus_->push_back(tmpTau);
            ACPFTauRef tmpRef(pfTaus_->back());
            pfTauMatching_->logConversion(reco::PFTauRef(pfTaus, index), tmpRef);
        }
    }

    /// set kinematic decays
    edm::Handle<SelectedKinematicDecayCollection> kinematicTaus;
    if (loadCollection(evt, kinematicTausTag_, kinematicTaus)) {
        *kinematicDecays_ = std::vector<ACFittedThreeProngDecay *>();
        *kinematicParticles_ = std::vector<ACFitParticle *>();
        for (SelectedKinematicDecayCollection::const_iterator decay = kinematicTaus->begin(); decay != kinematicTaus->end(); ++decay) {
            ACFittedThreeProngDecayConverter tmp(evt, *decay, kinematicParticleMatching_, kinematicParticles_, pfTauMatching_);
            ACFittedThreeProngDecay * tmpP = new ACFittedThreeProngDecay();
            *tmpP = tmp;
            kinematicDecays_->push_back(tmpP);
        }
    }
}

template <class T> void FinalTreeFiller::deleteVectorOfPointers(T * inVectorOfPointers) {
    if (!inVectorOfPointers) {
        return;
    }
    typename T::iterator i;
    for (i = inVectorOfPointers->begin(); i < inVectorOfPointers->end(); ++i) {
        if (!*i) continue;
        delete * i;
    }
}
template <class T> bool FinalTreeFiller::loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle) {
    iEvent.getByLabel(tag, handle);
    if (!handle.isValid()) {
        edm::LogError("FinalTreeFiller") << "No valid handle found for '" << tag << "'!";
        return false;
    }
    return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FinalTreeFiller);
