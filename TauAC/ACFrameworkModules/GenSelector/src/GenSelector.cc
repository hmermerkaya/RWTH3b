#include "../interface/GenSelector.h"

GenSelector::GenSelector(const edm::ParameterSet& iConfig):
candCollection_(iConfig.getParameter<edm::InputTag>("candCollection")),
decayType_(iConfig.getUntrackedParameter("decayType",std::string(""))),
motherPdgID_(iConfig.getUntrackedParameter("motherPdgID", 0))
{
    produces<int>("flag");//0=invalid, 1=valid

    produces<reco::GenParticleCollection>("genSignalDecay");//whole signal decay, fixed size and order: (0)mother, tau-, (unstable pi0's), stable daughters, tau+, (unstable pi0's), stable daughters, (quarks)
    produces<reco::GenParticleRefVector>("genSignalDecayRef");//store refs to input list to access signal particles in matchMap, same length/content as genSignalDecay, but originals not copied objects
    produces<reco::GenParticleRefVector>("genSignalQuarksRef");//only both VBF quarks, seperated for JetAnalysis. ref points to newly created genSignalDecay (not to initial collection)!!!
}


GenSelector::~GenSelector() {
}



// ------------ method called on each new Event  ------------
bool GenSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    cnt_++;

    genSignalPtr_ = std::auto_ptr<reco::GenParticleCollection>(new reco::GenParticleCollection); genSignalPtr_->reserve(13);
    genSignalRefPtr_ = std::auto_ptr<reco::GenParticleRefVector >(new reco::GenParticleRefVector);
    genSignalQuarksRefPtr_ = std::auto_ptr<reco::GenParticleRefVector >(new reco::GenParticleRefVector);

    reco::GenParticleCollection & genSignal = * genSignalPtr_;
    reco::GenParticleRefVector & genSignalRef = * genSignalRefPtr_;
    reco::GenParticleRefVector & genSignalQuarksRef = * genSignalQuarksRefPtr_;

    bool filterValue = false;
    if (iEvent.getByLabel(candCollection_, genCandidate_)) {//"genParticleCandidates"
        filterValue = storeEvent(iEvent, genSignal, genSignalRef);
    } else {
        LogTrace("GenSelector")<<"GenSelector::filter: Desired MC collection "<<candCollection_.label()<<" not found. Assuming to run on data and continue.";
        // place an empty collection into the event stream in case there is no generator information (as in data)
        // printf("no MC collection found\n");
    }

    edm::OrphanHandle<reco::GenParticleCollection> orphanSignal = iEvent.put(genSignalPtr_,"genSignalDecay");
    iEvent.put(genSignalRefPtr_,"genSignalDecayRef");
    genSignalQuarksRef = fillQuarkRefs(orphanSignal);
    iEvent.put(genSignalQuarksRefPtr_,"genSignalQuarksRef");


    std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
    int &flag = *flagPtr;
    if (filterValue) flag = 1;
    else flag = 0;
    iEvent.put(flagPtr,"flag");

    return filterValue;
}

// ------------ method called once each job just before starting event loop  ------------
void GenSelector::beginJob() {
    cnt_ = 0;
    cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void GenSelector::endJob() {
    float ratio = 0.0;
    if (cnt_!=0) ratio=(float)cntFound_/cnt_;
    edm::LogVerbatim("GenSelector")<<"--> [GenSelector] found at least one charged tau from "<<motherPdgID_<<". Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<ratio*100.0;
}

bool GenSelector::storeEvent(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef) {
    return storeEventAllTau(iEvent, collection, collectionRef, motherPdgID_);//store each tau followed by its stable daughters. if motherPdgID defined, only taus from this mother are kept. in case of VBF Higgs, also the quarks are stored.

    // place new generator selectors here
    edm::LogError("GenSelector")<<"GenSelector::storeEvent: The decay type "<<decayType_<<" (motherPdgID "<<motherPdgID_<<") is not implemented yet."<<std::endl;
    return false;
}
bool GenSelector::storeEventAllTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef, int motherPdgID) {
    bool gen = false;
    std::vector<reco::GenParticleRef> particlesRef;
    std::vector<reco::GenParticleRef> quarksRef;

    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index) {
        // store the given mother (if any)
        if (motherPdgID != 0) {
            if (abs(candIter->pdgId())==motherPdgID && candIter->numberOfDaughters()==3) {//i.e. Z->tautau(Z)
                if (particlesRef.size() != 0) {
                    edm::LogWarning("GenSelector")<<"GenSelector::storeEventAllTau:ERROR! More than one mother found. Skip it.";
                    continue;
                }                
                reco::GenParticleRef ref(genCandidate_, index);
                particlesRef.push_back(ref);
                //std::cout<<"GenSelector:: Found mother with pdgID"<<motherPdgID<<"."<<std::endl;
            }
        }
        // store quarks (e.g. from VBF).
        if (motherPdgID == 25) {
            if (abs(candIter->pdgId()) >= 1 && abs(candIter->pdgId()) <= 6 && candIter->numberOfDaughters()==1) {
                const reco::Candidate *genMother = (*candIter).mother();
                if (!genMother) continue;
                for (unsigned int k = 0; k!=genMother->numberOfDaughters(); k++) {
                    if (abs(genMother->daughter(k)->pdgId())==25 && genMother->numberOfDaughters()==3) {
                        reco::GenParticleRef ref(genCandidate_, index);
                        quarksRef.push_back(ref);
                    }
                }
            }
        }
        // store taus
        if (abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()>1) {//ignore tau->tau
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            genMother = ((*candIter).mother())->mother();
            if (!genMother) continue;
            if (motherPdgID != 0) {
                if (particlesRef.size() < 1) {
                    LogTrace("GenSelector")<<"GenSelector::storeEventAllTau:ERROR! Mother was not found before occurance of a tau. Skip this tau.";
                    continue;
                }
                if (genMother != particlesRef.at(0).get()) {
                    LogTrace("GenSelector")<<"GenSelector::storeEventAllTau: Skipping tau not from "<<motherPdgID<<" but from "<<genMother->pdgId()<<".";
                    continue;
                }
            }
            reco::GenParticleRef tauRef(genCandidate_, index);

            std::vector<reco::GenParticleRef> daughtersRef;
            findDescendents(tauRef, daughtersRef, 2, 111);// add unstable pi0's
            findDescendents(tauRef, daughtersRef, 1);//add stable daughters
            //std::cout<<"GenSelector:: found tau from "<<genMother->pdgId()<<" decaying into "<<daughtersRef.size()<<" daughters."<<std::endl;
            particlesRef.push_back(tauRef);
            particlesRef.insert(particlesRef.end(), daughtersRef.begin(), daughtersRef.end());
        }
    }
    if (particlesRef.size()>0) {
        // printf("GenSelector:: found %lu particles (either taus or their stable descendents).\n", particlesRef.size());
        cntFound_++;
        gen = true;

        for (std::vector<reco::GenParticleRef>::const_iterator idr = particlesRef.begin(); idr!= particlesRef.end(); ++idr) {
            collection.push_back(**idr);
            collectionRef.push_back(*idr);
        }
        if (quarksRef.size()>0) {        
            sort(quarksRef.begin(),quarksRef.end(), ACTools::cmpPz<const reco::GenParticleRef>);
            for (std::vector<reco::GenParticleRef>::const_iterator idr = quarksRef.begin(); idr!= quarksRef.end(); ++idr) {
                collection.push_back(**idr);
                collectionRef.push_back(*idr);
            }
        }
    } else {
        // std::cout<<"evt"<<iEvent.id().event()<<" GenSelector::storeEventAllTau: no tau decay found in event type: "<<decayType_<<std::endl;
        gen = false;
    }

    return gen;
}
void GenSelector::findDescendents(const reco::GenParticleRef& base, std::vector<reco::GenParticleRef> & descendents, int status, int pdgId) {
    const reco::GenParticleRefVector& daughterRefs = base->daughterRefVector();
    for (reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
        if ((*idr)->status() == status && (!pdgId || abs((*idr)->pdgId()) == pdgId)) descendents.push_back(*idr);
        else findDescendents(*idr, descendents, status, pdgId);
    }
}
reco::GenParticleRefVector GenSelector::fillQuarkRefs(const edm::OrphanHandle<reco::GenParticleCollection> & orphanSignal) {
    reco::GenParticleRefVector refs;
    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator iter = orphanSignal->begin(); iter != orphanSignal->end(); ++iter, ++index) {
        reco::GenParticleRef ref(orphanSignal, index);
        refs.push_back(ref);
    }

    return refs;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSelector);
