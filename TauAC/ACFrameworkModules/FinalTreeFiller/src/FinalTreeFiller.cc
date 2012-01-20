#include "ACFrameworkModules/FinalTreeFiller/interface/FinalTreeFiller.h"

FinalTreeFiller::FinalTreeFiller(const edm::ParameterSet& iConfig):
genSignalTag_(iConfig.getParameter<edm::InputTag>("genSignal")),
genSignalRefTag_(iConfig.getParameter<edm::InputTag>("genSignalRef")),
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
pfTauDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("pfTauDiscriminators")),
flags_(iConfig.getParameter< std::vector<std::string> >("flags")),
pileUpDistributionFileMC_(iConfig.getUntrackedParameter<std::string >("pileUpDistributionFileMC", "")),
pileUpDistributionHistMC_(iConfig.getUntrackedParameter<std::string >("pileUpDistributionHistMC", "")),
pileUpDistributionFileData_(iConfig.getUntrackedParameter<std::string >("pileUpDistributionFileData", "")),
pileUpDistributionHistData_(iConfig.getUntrackedParameter<std::string >("pileUpDistributionHistData", ""))
{
    hltChanged_ = true;
    eventTree_ = 0;
    eventInfo_ = 0;
    eventWeight_ = 0;
    eventGlobals_ = 0;
    trigger_ = 0;
    offlinePV_ = 0;
    reducedPV_ = 0;
    muons_ = 0;
    electrons_ = 0;
    fittedThreeProngParticles_ = 0;
    tauDecays_ = 0;
    generator_ = 0;
    genTauDecays_ = 0;
    pfJets_ = 0;
    pfTaus_ = 0;
    pileup_ = 0;
    lumiWeights_ = 0;
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
    deleteVectorOfPointers(fittedThreeProngParticles_);
    deleteVectorOfPointers(tauDecays_);
    deleteVectorOfPointers(generator_);
    deleteVectorOfPointers(genTauDecays_);
    deleteVectorOfPointers(pfJets_);
    deleteVectorOfPointers(pfTaus_);
    deleteVectorOfPointers(pileup_);
}
void FinalTreeFiller::beginJob() {
    evtCnt_ = 0;
    runCnt_ = 0;

    kinematicParticleMatching_ = new MCTruthMatching(chargedTauDaughterMatchMapTag_);
    conversionLogPFTau_ = new PFTauConversionLog();
    conversionLogPFJet_ = new PFJetConversionLog();

    eventTree_ = fileService_->make<TTree>("TauACEvent", "TauAC tree filled each event");
    //runTree_ = fileService_->make<TTree>("TauACRun", "TauAC tree filled each run");

    eventInfo_ = new ACEventInfo();
    eventTree_->Branch("ACEventInfo", &eventInfo_, 32000, 0);

    eventWeight_ = new ACEventWeight();
    eventTree_->Branch("ACEventWeight", &eventWeight_, 32000, 0);
    
    eventGlobals_ = new ACEventGlobals();
    eventTree_->Branch("ACEventGlobals", &eventGlobals_, 32000, 0);

    trigger_ = new ACTrigger();
    eventTree_->Branch("ACTrigger", &trigger_, 32000, 0);

    offlinePV_ = new std::vector<ACVertex *>();
    eventTree_->Branch("offlinePV", &offlinePV_, 32000, 0);
    reducedPV_ = new std::vector<ACVertex *>();
    eventTree_->Branch("reducedPV", &reducedPV_, 32000, 0);

    generator_ = new std::vector<ACGenParticle *>();
    eventTree_->Branch("generator", &generator_, 32000, 0);

    genTauDecays_ = new std::vector<ACGenDecay *>();
    eventTree_->Branch("genTauDecays", &genTauDecays_, 32000, 0);

    muons_ = new std::vector<ACParticle *>();
    eventTree_->Branch("Muons", &muons_, 32000, 0);

    electrons_ = new std::vector<ACParticle *>();
    eventTree_->Branch("Electrons", &electrons_, 32000, 0);

    fittedThreeProngParticles_ = new std::vector<ACFitParticle *>();
    eventTree_->Branch("ACFittedThreeProngParticles", &fittedThreeProngParticles_, 32000, 0);

    tauDecays_ = new std::vector<ACFittedThreeProngDecay *>();
    eventTree_->Branch("ACFittedThreeProngDecays", &tauDecays_, 32000, 0);

    pfJets_ = new std::vector<ACJet *>();
    eventTree_->Branch("PFJets", &pfJets_, 32000, 0);

    pfTaus_ = new std::vector<ACPFTau *>();
    eventTree_->Branch("PFTaus", &pfTaus_, 32000, 0);
    
    pileup_ = new std::vector<ACPileupInfo *>();
    eventTree_->Branch("ACPileupInfo", &pileup_, 32000, 0);
    
    /// initialize pileup reweighting (it might be useful to renew the initialization for every run!?)
    if (pileUpDistributionFileMC_ != "" && pileUpDistributionFileData_ != "") {
        lumiWeights_ = new edm::Lumi3DReWeighting(pileUpDistributionFileMC_, pileUpDistributionFileData_, pileUpDistributionHistMC_, pileUpDistributionHistData_, "");
        lumiWeights_->weight3D_init(1.0);
    }
}
void FinalTreeFiller::endJob() {
    kinematicParticleMatching_->printOutro();

    /// clear memory
    delete kinematicParticleMatching_;
    delete conversionLogPFTau_;
    delete conversionLogPFJet_;
    delete eventInfo_;
    delete eventWeight_;
    delete eventGlobals_;
    delete trigger_;
    delete offlinePV_;
    delete reducedPV_;
    delete muons_;
    delete electrons_;
    delete fittedThreeProngParticles_;
    delete tauDecays_;
    delete generator_;
    delete genTauDecays_;
    delete pfJets_;
    delete pfTaus_;
    delete pileup_;
    delete lumiWeights_;
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
    conversionLogPFTau_->nextEvt();
    conversionLogPFJet_->nextEvt();
    
    /// set event information
    edm::Handle<reco::GenParticleCollection> genStatus;
    /// check whether generator information is present in the event
    if (evt.getByLabel(edm::InputTag("genParticles"), genStatus)) {
        *eventInfo_ = ACEventInfo(evt.id().event(), evt.id().luminosityBlock(), evt.id().run(), "MC");
    } else {
        *eventInfo_ = ACEventInfo(evt.id().event(), evt.id().luminosityBlock(), evt.id().run(), "data");
    }
    for (std::vector<std::string>::iterator flag=flags_.begin(); flag!=flags_.end(); ++flag) {
		edm::Handle<int> flagHandle;
        if (loadCollection(evt, edm::InputTag(*flag, "flag"), flagHandle)) {
            eventInfo_->storeEDMFilterResult(*flag, *flagHandle);
        }
	}
    
    /// calculate and store event weights
    double eventWeight = 1.0;
    if (eventInfo_->type() == "MC") {
        const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&evt);
        eventWeight = lumiWeights_->weight3D(*iEventB);
    }
    *eventWeight_ = ACEventWeight(eventWeight);
    
    /// set event globals
    edm::Handle<reco::PFMETCollection> pfMET;
    edm::Handle<reco::METCollection> tcMET;
    if (loadCollection(evt, pfMETTag_, pfMET) && loadCollection(evt, tcMETTag_, tcMET)) {
        *eventGlobals_ = ACEventGlobals(TVector3(pfMET->begin()->px(), pfMET->begin()->py(), 0.0), TVector3(tcMET->begin()->px(), tcMET->begin()->py(), 0.0), pfMET->begin()->sumEt(), tcMET->begin()->sumEt());
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
            *trigger_ = ACTrigger(HLTCP_.tableName(), hltPathMap);
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
    if (loadCollection(evt, genSignalTag_, genCands, true)) {
        *generator_ = std::vector<ACGenParticle *>();
        *genTauDecays_ = std::vector<ACGenDecay *>();
        ACGenDecayRef latestGenDecayRef;
        unsigned int index = 0;
        for (reco::GenParticleCollection::const_iterator candidate = genCands->begin(); candidate != genCands->end(); ++candidate, index++) {
            ACGenParticleConverter tmp(*candidate);
            ACGenParticle * tmpP = new ACGenParticle();
            *tmpP = tmp;
            generator_->push_back(tmpP);
            ACGenParticleRef tmpRef(generator_->back());
            kinematicParticleMatching_->logConversion(reco::GenParticleRef(genCands, index), tmpRef);

            // store mother/daughter relation. additional loop needed!!!
            if (candidate->motherRefVector().size()==1) {
                const reco::GenParticleRef & genMother = candidate->motherRef(0);//WARNING!!! mother is stored in the initial collection. not the one I had to use for the matching!!!
                // ugly: need the references to the initial collection here, due to mother() function
                // dangerous: or should we modify the mother link in the selected collection?! is this possible at all? yes, using CompositeRefCandidate::addDaughter() and co
                edm::Handle<reco::GenParticleRefVector> genCandRefs;
                if (loadCollection(evt, genSignalRefTag_, genCandRefs, true)) {
                    unsigned int refIndex = 0;
                    for (reco::GenParticleRefVector::const_iterator candidateRef = genCandRefs->begin(); candidateRef != genCandRefs->end(); ++candidateRef, refIndex++) {
                        if (genMother == *candidateRef){//as not all unstable particles are stored, a mother is not always found!
                            //printf("for particle with pdgID %i found mother with pdgID %i at refIndex %i. size is %lu!\n", tmpRef->pdgId(), genMother->pdgId(), refIndex, generator_->size());
                            if (refIndex >= generator_->size()) {
                                printf("FinalTreeFiller::storeEvent:ERROR! Bad order in generator collection. Mother stored after daughter.\n");
                                throw 404;
                            }    
                            ACGenParticleRef motherRef(generator_->at(refIndex));
                            ACGenParticleConverter * tmpConverter = static_cast<ACGenParticleConverter*>(generator_->back());
                            tmpConverter->setMother(motherRef);
                            tmpConverter = static_cast<ACGenParticleConverter*>(generator_->at(refIndex));
                            tmpConverter->addDaughter(tmpRef);
                            break;
                        }
                    }
                }
            }

            // store every occuring tau decay, expect the tau to be followed by its daughters
            if (abs(candidate->pdgId())==15) {
                ACGenDecayConverter tmpDecay(tmpRef);
                ACGenDecay * tmpDecayP = new ACGenDecay();
                *tmpDecayP = tmpDecay;
                genTauDecays_->push_back(tmpDecayP);
            } else {
                // skip all particles before the first tau.
                if (genTauDecays_->size()<1) continue;
                // skip quarks
                if (abs(candidate->pdgId()) >= 1 && abs(candidate->pdgId()) <= 6) continue;
                ACGenDecayConverter * tmpConverter = static_cast<ACGenDecayConverter*>(genTauDecays_->back());
                if (tmpConverter) {
                    tmpConverter->addDaughter(tmpRef);
                }
            }
            
            // store the reference to the generator tau decay that contains the current particle
            ACGenParticleConverter * tmpConverter = static_cast<ACGenParticleConverter*>(generator_->back());
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
        unsigned int index = 0;
        for (reco::PFJetCollection::const_iterator ijet = pfJets->begin(); ijet != pfJets->end(); ++ijet, index++) {
            ACJetConverter tmp(*ijet);
            ACJet * tmpJet = new ACJet();
            *tmpJet = tmp;
            pfJets_->push_back(tmpJet);
            ACJetRef tmpRef(pfJets_->back());
            conversionLogPFJet_->logConversion(reco::PFJetRef(pfJets, index), tmpRef);
        }
    }
    
    /// set pf tau collection
    edm::Handle<reco::PFTauCollection> pfTaus;
    if (loadCollection(evt, pfTauTag_, pfTaus)) {
        *pfTaus_ = std::vector<ACPFTau *>();
        unsigned int index = 0;
        for (reco::PFTauCollection::const_iterator itau = pfTaus->begin(); itau != pfTaus->end(); ++itau, index++) {
            ACPFTauConverter tmp(evt, reco::PFTauRef(pfTaus, index), pfTauDiscriminatorTags_, conversionLogPFJet_);
            ACPFTau * tmpTau = new ACPFTau();
            *tmpTau = tmp;
            pfTaus_->push_back(tmpTau);
            ACPFTauRef tmpRef(pfTaus_->back());
            conversionLogPFTau_->logConversion(reco::PFTauRef(pfTaus, index), tmpRef);
        }
    }
    
    /// set kinematic decays
    edm::Handle<SelectedKinematicDecayCollection> kinematicTaus;
    if (loadCollection(evt, kinematicTausTag_, kinematicTaus)) {
        *tauDecays_ = std::vector<ACFittedThreeProngDecay *>();
        *fittedThreeProngParticles_ = std::vector<ACFitParticle *>();
        for (SelectedKinematicDecayCollection::const_iterator decay = kinematicTaus->begin(); decay != kinematicTaus->end(); ++decay) {
            ACFittedThreeProngDecayConverter tmp(evt, *decay, kinematicParticleMatching_, fittedThreeProngParticles_, conversionLogPFTau_);
            ACFittedThreeProngDecay * tmpP = new ACFittedThreeProngDecay();
            *tmpP = tmp;
            tauDecays_->push_back(tmpP);
        }
    }
    
    /// set pileup information
    edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;
    if (loadCollection(evt, pileupInfoTag_, pileupInfo, true)) {
        *pileup_ = std::vector<ACPileupInfo *>();
        /// one PileupSummaryInfo for each of the beam crossings (allows to retrieve information about the out-of-time pileup)
        for (std::vector<PileupSummaryInfo>::const_iterator PVI = pileupInfo->begin(); PVI != pileupInfo->end(); ++PVI) {
            ACPileupConverter tmp2(*PVI);
            ACPileupInfo * tmpPileup = new ACPileupInfo();
            *tmpPileup = tmp2;
            pileup_->push_back(tmpPileup);
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(FinalTreeFiller);
