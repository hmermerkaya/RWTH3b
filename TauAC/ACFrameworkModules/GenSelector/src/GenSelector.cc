#include "../interface/GenSelector.h"

GenSelector::GenSelector(const edm::ParameterSet& iConfig):
candCollection_(iConfig.getParameter<edm::InputTag>("candCollection")),
decayType_(iConfig.getUntrackedParameter("decayType",std::string("VBFH"))),
motherPdgID_(iConfig.getUntrackedParameter("motherPdgID", 0))
{
    produces<int>("flag");//0=invalid, 1=valid

    produces<reco::GenParticleCollection>("genSignalDecay");//whole signal decay, fixed size and order: (0)higgs, (1)tauM, (2..4)piM, (5)nuM, (6)tauP, (7..9)piP, (10)nuP, (11..12)quarks = 13 particles
    produces<reco::GenParticleRefVector>("genSignalDecayRef");//store refs to input list to access signal particles in matchMap
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
        LogTrace("MCCollectionNotFound")<<"GenSelector::filter: Desired MC collection "<<candCollection_.label()<<" not found. Assuming to run on data and continue.";
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
    if (decayType_=="TauPair") printf("--> [GenSelector] found at least two opp. charged taus from %d. Efficiency: %d/%d = %.2f%%\n", motherPdgID_, cntFound_, cnt_, ratio*100.0);
    else if (decayType_!="unknown") printf("--> [GenSelector] found at least two opp. charged taus from %s. Efficiency: %d/%d = %.2f%%\n", decayType_.c_str(), cntFound_, cnt_, ratio*100.0);
    else printf("--> [GenSelector] found at least one unknown tau decay. Efficiency: %d/%d = %.2f%%\n", cntFound_, cnt_, ratio*100.0);
}

bool GenSelector::storeEvent(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef) {
    if (decayType_=="VBFH3pr") return storeEventVBFH(iEvent, collection, collectionRef);
    if (decayType_=="W3pr") return storeEventSingleTau(iEvent, collection, collectionRef);
    if (decayType_=="TauPair" && motherPdgID_!=0) return storeTauPairEvt(iEvent, collection, collectionRef, motherPdgID_);
    if (decayType_=="unknown") return storeEventAllTau(iEvent, collection, collectionRef);//search for all taus

    std::cout<<"The decay type "<<decayType_<<" (motherPdgID "<<motherPdgID_<<") is not implemented yet."<<std::endl;
    return false;
}
bool GenSelector::storeTauPairEvt(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef, const int & motherPdgID) {
    bool gen = false;
    std::vector<const reco::GenParticle* > mother; // Z or H, W, ...
    std::vector<const reco::GenParticle* > tauP; // charge+
    std::vector<const reco::GenParticle* > tauM; // charge-
    std::vector<const reco::GenParticle* > daughtersP; // all stable daughters
    std::vector<const reco::GenParticle* > daughtersM; // all stable daughters

    std::vector<reco::GenParticleRef> motherRef;
    std::vector<reco::GenParticleRef> tauPRef; // charge+
    std::vector<reco::GenParticleRef> tauMRef; // charge-
    std::vector<reco::GenParticleRef> daughtersPRef;
    std::vector<reco::GenParticleRef> daughtersMRef;

    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index) {
        if (abs(candIter->pdgId())==motherPdgID && candIter->numberOfDaughters()==3) {//i.e. Z->tautau(Z)
            mother.push_back(&*candIter);
            reco::GenParticleRef ref(genCandidate_, index);
            motherRef.push_back(ref);
            // std::cout<<"GenSelector:: Found mother with pdgID"<<motherPdgID<<"."<<std::endl;
        }
        if (abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()>1) {//ignore tau->tau
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            genMother = ((*candIter).mother())->mother();
            if (!genMother && mother.size()>0) continue;
            if (genMother == mother.at(0)) {
                // std::cout<<"GenSelector:: found tau from "<<motherPdgID<<"."<<std::endl;
                if (candIter->pdgId()==15) {
                    tauM.push_back(&*candIter);
                    reco::GenParticleRef ref(genCandidate_, index);
                    tauMRef.push_back(ref);
                    findDescendents(ref, daughtersMRef, 1);
                } else {
                    tauP.push_back(&*candIter);
                    reco::GenParticleRef ref(genCandidate_, index);
                    tauPRef.push_back(ref);
                    findDescendents(ref, daughtersPRef, 1);
                }
            } // else std::cout<<"GenSelector:: found tau not from "<<motherPdgID<<"."<<std::endl;
        }
    }
    // printf("GenSelector:: found %lu + %lu stable tau descendents.\n", daughtersMRef.size(), daughtersPRef.size());
    if (mother.size()!=1) {
        if (mother.size()<1) std::cout<<"No mother with pdgID "<<motherPdgID<<" detected.";
        if (mother.size()>1) std::cout<<"More than one mother with pdgID "<<motherPdgID<<" detected.";
        std::cout<<" Skip event."<<std::endl;
        return false;
    }
    if (tauP.size()==1 && tauM.size()==1 && daughtersPRef.size()>1 && daughtersMRef.size()>1) {
        cntFound_++;
        gen = true;

        // fixed order
        collection.push_back(*(mother.at(0)));
        collection.push_back(*(tauM.at(0)));
        for (std::vector<reco::GenParticleRef>::const_iterator idr = daughtersMRef.begin(); idr!= daughtersMRef.end(); ++idr) daughtersM.push_back(&**idr);
        for (std::vector<const reco::GenParticle* >::const_iterator idr = daughtersM.begin(); idr!= daughtersM.end(); ++idr) collection.push_back(**idr);
        collection.push_back(*(tauP.at(0)));
        for (std::vector<reco::GenParticleRef>::const_iterator idr = daughtersPRef.begin(); idr!= daughtersPRef.end(); ++idr) daughtersP.push_back(&**idr);
        for (std::vector<const reco::GenParticle* >::const_iterator idr = daughtersP.begin(); idr!= daughtersP.end(); ++idr) collection.push_back(**idr);

        collectionRef.push_back(motherRef.at(0));
        collectionRef.push_back(tauMRef.at(0));
        for (std::vector<reco::GenParticleRef>::const_iterator idr = daughtersMRef.begin(); idr!= daughtersMRef.end(); ++idr) collectionRef.push_back(*idr);
        collectionRef.push_back(tauPRef.at(0));
        for (std::vector<reco::GenParticleRef>::const_iterator idr = daughtersPRef.begin(); idr!= daughtersPRef.end(); ++idr) collectionRef.push_back(*idr);
    } else {
        std::cout<<"evt"<<iEvent.id().event()<<" GenSelector::storeEventZSMBr: "<<motherPdgID<<"->tau+tau->SMBr not found\t";
        //        std::cout<<"#of tau "<<tauCand->size()<<"#of nu "<<nuCand->size()<<"#of pi "<<piCand->size()<<std::endl;
        //    std::cout<<"evtType"<<std::endl;
        std::cout<<"sizes: taus "<<tauM.size()<<"+"<<tauP.size()<<", daughters "<<daughtersM.size()<<"+"<<daughtersP.size()<<std::endl;
        std::cout<<"gen-";
        gen = false;
    }

    return gen;
}
bool GenSelector::storeEventSingleTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef) {
    bool gen = false;
    std::vector<const reco::GenParticle* > mother;//W
    std::vector<const reco::GenParticle* > tau;
    std::vector<const reco::GenParticle* > pi;
    std::vector<const reco::GenParticle* > nu;

    std::vector<reco::GenParticleRef> motherRef;
    std::vector<reco::GenParticleRef> tauRef;
    std::vector<reco::GenParticleRef> piRef;
    std::vector<reco::GenParticleRef> nuRef;

    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index) {
        if (abs(candIter->pdgId())==24 && candIter->numberOfDaughters()==3) {//W->tauNu(W)
            mother.push_back(&*candIter);
            reco::GenParticleRef ref(genCandidate_, index);
            motherRef.push_back(ref);
        }
        if (abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()==4) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            genMother = ((*candIter).mother())->mother();
            if (!genMother && mother.size()>0) continue;
            if (genMother == mother.at(0)) {
                int cntDaughterChecked = 0;
                for (unsigned int k = 0; k!=candIter->numberOfDaughters(); k++) {
                    if (abs(candIter->daughter(k)->pdgId())!=16 && abs(candIter->daughter(k)->pdgId())!=211) break;
                    cntDaughterChecked++;
                }
                if (cntDaughterChecked!=4) continue;
                tau.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                tauRef.push_back(ref);
            } // else std::cout<<"GenSelector:: found other tau decaying into 3pi+nu."<<std::endl;
        }
        // std::cout<<"taus: "<<tauM.size()+tauP.size()<<std::endl;
        if (tau.size()==0) continue;
        if (abs(candIter->pdgId()) == 16) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            if (candIter->mother() == tau.at(0)) {
                nu.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                nuRef.push_back(ref);
            }
        }
        if (abs(candIter->pdgId()) == 211) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            if (candIter->mother() == tau.at(0)) {
                pi.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                piRef.push_back(ref);
            }
        }
    }
    if (mother.size()!=1) {
        if (mother.size()<1) std::cout<<"No W detected.";
        if (mother.size()>1) std::cout<<"More than one W detected.";
        std::cout<<" Skip event."<<std::endl;
        return false;
    }
    if (tau.size()==1 && nu.size()==1 && pi.size()==3) {
        cntFound_++;
        gen = true;

        // fixed order and size!
        collection.push_back(*(mother.at(0)));
        collection.push_back(*(tau.at(0)));
        sort(pi.begin(),pi.end(), ACTools::cmpPt<const reco::GenParticle*>);
        for (unsigned int i=0; i!=pi.size(); i++) collection.push_back(*(pi.at(i)));
        collection.push_back(*(nu.at(0)));

        collectionRef.push_back(motherRef.at(0));
        collectionRef.push_back(tauRef.at(0));
        sort(piRef.begin(),piRef.end(), ACTools::cmpPt<reco::GenParticleRef>);
        for (unsigned int i=0; i!=piRef.size(); i++) collectionRef.push_back(piRef.at(i));
        collectionRef.push_back(nuRef.at(0));
        // printf("nuM = (%8.6f, %8.6f, %8.6f)\n", nuMRef.at(0)->vx(), nuMRef.at(0)->vy(), nuMRef.at(0)->vz());
        // printf("nuP = (%8.6f, %8.6f, %8.6f)\n", nuPRef.at(0)->vx(), nuPRef.at(0)->vy(), nuPRef.at(0)->vz());
    } else {
        std::cout<<"evt"<<iEvent.id().event()<<" GenSelector::storeEventSingleTau: W->tau+nu->3prong not found\t";
        // std::cout<<"#of tau "<<tauCand->size()<<"#of nu "<<nuCand->size()<<"#of pi "<<piCand->size()<<std::endl;
        // std::cout<<"evtType"<<std::endl;
        std::cout<<"sizes: nus "<<nu.size()<<", pis "<<pi.size()<<std::endl;
        std::cout<<"gen-";
        gen = false;
    }

    return gen;
}
bool GenSelector::storeEventAllTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef) {
    bool gen = false;
    std::vector<reco::GenParticleRef> particlesRef;

    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index) {
        if (abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()>1) {//ignore tau->tau
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            genMother = ((*candIter).mother())->mother();
            if (!genMother) continue;
            reco::GenParticleRef tauRef(genCandidate_, index);

            std::vector<reco::GenParticleRef> daughtersRef;
            findDescendents(tauRef, daughtersRef, 1);
            std::cout<<"GenSelector:: found tau from "<<genMother->pdgId()<<" decaying into "<<daughtersRef.size()<<" daughters."<<std::endl;
            particlesRef.push_back(tauRef);
            particlesRef.insert(particlesRef.end(), daughtersRef.begin(), daughtersRef.end());
        }
    }
    if (particlesRef.size()>0) {
        // printf("GenSelector:: found %lu particles (either taus or their stable descendents).\n", particlesRef.size());
        cntFound_++;
        gen = true;

        // WARNING: no mother at the beginning!
        for (std::vector<reco::GenParticleRef>::const_iterator idr = particlesRef.begin(); idr!= particlesRef.end(); ++idr) {
            collection.push_back(**idr);
            collectionRef.push_back(*idr);
        }
    } else {
        // std::cout<<"evt"<<iEvent.id().event()<<" GenSelector::storeEventAllTau: no tau decay found in event type: "<<decayType_<<std::endl;
        gen = false;
    }

    return gen;
}
bool GenSelector::storeEventVBFH(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef) {
    bool gen = false;
    std::vector<const reco::GenParticle* > mother; // higgs
    std::vector<const reco::GenParticle* > tauP; // charge+
    std::vector<const reco::GenParticle* > tauM; // charge-
    std::vector<const reco::GenParticle* > piP;
    std::vector<const reco::GenParticle* > piM;
    std::vector<const reco::GenParticle* > nuP;
    std::vector<const reco::GenParticle* > nuM;
    std::vector<const reco::GenParticle* > quarks;

    std::vector<reco::GenParticleRef> motherRef;
    std::vector<reco::GenParticleRef> tauPRef; // charge+
    std::vector<reco::GenParticleRef> tauMRef; // charge-
    std::vector<reco::GenParticleRef> piPRef;
    std::vector<reco::GenParticleRef> piMRef;
    std::vector<reco::GenParticleRef> nuPRef;
    std::vector<reco::GenParticleRef> nuMRef;
    std::vector<reco::GenParticleRef> quarksRef;


    // for (reco::CandidateView::const_iterator candIter = candView->begin(); candIter != candView->end(); ++candIter) {
    unsigned index = 0;
    for (reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index) {
        if (abs(candIter->pdgId()) >= 1 && abs(candIter->pdgId()) <= 6 && candIter->numberOfDaughters()==1) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            for (unsigned int k = 0; k!=genMother->numberOfDaughters(); k++) {
                if (abs(genMother->daughter(k)->pdgId())==25 && genMother->numberOfDaughters()==3) {
                    quarks.push_back(&*candIter);
                    reco::GenParticleRef ref(genCandidate_, index);
                    quarksRef.push_back(ref);
                }
            }
        }
        // std::cout<<"quarks: "<<quarks.size()<<std::endl;
        if (abs(candIter->pdgId())==25 && candIter->numberOfDaughters()==3) {//h0->2tau+h0
            mother.push_back(&*candIter);
            // edm::Ref<C> myRef(const edm::Handle<C> & collection, size_type indexIntoCollection)
            reco::GenParticleRef ref(genCandidate_, index);
            motherRef.push_back(ref);
        }
        // std::cout<<"higgs: "<<mother.size()<<std::endl;
        if (abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()==4) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            genMother = ((*candIter).mother())->mother();
            if (!genMother && mother.size()>0) continue;
            if (genMother == mother.at(0)) {
                int cntDaughterChecked = 0;
                for (unsigned int k = 0; k!=candIter->numberOfDaughters(); k++) {
                    if (abs(candIter->daughter(k)->pdgId())!=16 && abs(candIter->daughter(k)->pdgId())!=211) break;
                    cntDaughterChecked++;
                }
                if (cntDaughterChecked!=4) continue;
                if (candIter->pdgId()==15) {
                    tauM.push_back(&*candIter);
                    reco::GenParticleRef ref(genCandidate_, index);
                    tauMRef.push_back(ref);
                } else {
                    tauP.push_back(&*candIter);
                    reco::GenParticleRef ref(genCandidate_, index);
                    tauPRef.push_back(ref);
                }
            } // else std::cout<<"GenSelector:: found other tau decaying into 3pi+nu."<<std::endl;
        }
        // std::cout<<"taus: "<<tauM.size()+tauP.size()<<std::endl;
        if (tauM.size()==0 && tauP.size()==0) continue;
        if (abs(candIter->pdgId()) == 16) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            if (candIter->mother()->pdgId() == 15) if (candIter->mother() == tauM.at(0)) {
                nuM.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                nuMRef.push_back(ref);
            }
            if (candIter->mother()->pdgId() ==-15) if (candIter->mother() == tauP.at(0)) {
                nuP.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                nuPRef.push_back(ref);
            }
        }
        if (abs(candIter->pdgId()) == 211) {
            const reco::Candidate *genMother = (*candIter).mother();
            if (!genMother) continue;
            if (candIter->mother()->pdgId() == 15) if (candIter->mother() == tauM.at(0)) {
                piM.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                piMRef.push_back(ref);
            }
            if (candIter->mother()->pdgId() ==-15) if (candIter->mother() == tauP.at(0)) {
                piP.push_back(&*candIter);
                reco::GenParticleRef ref(genCandidate_, index);
                piPRef.push_back(ref);
            }
        }
    }
    if (mother.size()!=1) {
        if (mother.size()<1) std::cout<<"No higgs detected.";
        if (mother.size()>1) std::cout<<"More than one higgs detected.";
        std::cout<<" Skip event."<<std::endl;
        return false;
    }
    //    for (unsigned int i=0; i!=mother.size(); i++) {
    //        for (unsigned int j = 0; j!=mother.at(i)->numberOfDaughters()-1; j++) {//skip repeated higgs
    //            if (mother.at(i)->daughter(j)->numberOfDaughters()!=1) {
    //                std::cout<<"ERROR: tau not listed twice. Skip it."<<std::endl;
    //                continue;
    //            }
    //            if (abs(mother.at(i)->daughter(j)->daughter(0)->pdgId()) == 15 && mother.at(i)->daughter(j)->daughter(0)->numberOfDaughters()==4) {
    //                int cntDaughterChecked = 0;
    //                int charge = 0;
    //                if (mother.at(i)->daughter(j)->daughter(0)->pdgId()==15) charge = -1;
    //                else charge = 1;
    //                for (unsigned int k = 0; k!=mother.at(i)->daughter(j)->daughter(0)->numberOfDaughters(); k++) {
    //                    if (abs(mother.at(i)->daughter(j)->daughter(0)->daughter(k)->pdgId())==16) {
    //                        if (charge==-1) nuM.push_back(mother.at(i)->daughter(j)->daughter(0)->daughter(k));
    //                        else nuP.push_back(mother.at(i)->daughter(j)->daughter(0)->daughter(k));
    //                    }
    //                    if (abs(mother.at(i)->daughter(j)->daughter(0)->daughter(k)->pdgId())==211) {
    //                        if (charge==-1) piM.push_back(mother.at(i)->daughter(j)->daughter(0)->daughter(k));
    //                        else piP.push_back(mother.at(i)->daughter(j)->daughter(0)->daughter(k));
    //                    }
    //                    if (abs(mother.at(i)->daughter(j)->daughter(0)->daughter(k)->pdgId())!=16 && abs(mother.at(i)->daughter(j)->daughter(0)->daughter(k)->pdgId())!=211) {
    //                        break;
    //                    }
    //                    cntDaughterChecked++;
    //                }
    //                //                            std::cout<<"cntDaughterCheck: "<<cntDaughterChecked<<std::endl;
    //                if (cntDaughterChecked==4) {
    //                    if (charge==-1) tauM.push_back(mother.at(i)->daughter(j)->daughter(0));
    //                    else tauP.push_back(mother.at(i)->daughter(j)->daughter(0));
    //                }
    //            }
    //        }
    //        for (unsigned int j=0; j!=mother.at(i)->numberOfMothers(); j++) {
    //            if (abs(mother.at(i)->mother(j)->pdgId())>=1 && abs(mother.at(i)->mother(j)->pdgId())<=6) {
    //                quarks.push_back(mother.at(i)->mother(j));
    //            } else {
    //                std::cout<<"Bad higgs mother!"<<std::endl;
    //                continue;
    //            }
    //        }
    //    }
    if (tauP.size()==1 && tauM.size()==1 && nuP.size()==1 && nuM.size()==1 && piP.size()==3 && piM.size()==3 && quarks.size()==2) {
        // std::cout<<"H->tau+tau->2x3prong found";
        cntFound_++;
        gen = true;

        // fixed order and size!
        collection.push_back(*(mother.at(0)));
        collection.push_back(*(tauM.at(0)));
        sort(piM.begin(),piM.end(), ACTools::cmpPt<const reco::GenParticle*>);
        for (unsigned int i=0; i!=piM.size(); i++) collection.push_back(*(piM.at(i)));
        collection.push_back(*(nuM.at(0)));
        collection.push_back(*(tauP.at(0)));
        sort(piP.begin(),piP.end(), ACTools::cmpPt<const reco::GenParticle*>);
        for (unsigned int i=0; i!=piP.size(); i++) collection.push_back(*(piP.at(i)));
        collection.push_back(*(nuP.at(0)));
        sort(quarks.begin(),quarks.end(), ACTools::cmpPz<const reco::GenParticle*>);
        for (unsigned int i=0; i!=quarks.size(); i++) collection.push_back(*(quarks.at(i)));


        collectionRef.push_back(motherRef.at(0));
        collectionRef.push_back(tauMRef.at(0));
        sort(piMRef.begin(),piMRef.end(), ACTools::cmpPt<reco::GenParticleRef>);
        for (unsigned int i=0; i!=piMRef.size(); i++) collectionRef.push_back(piMRef.at(i));
        collectionRef.push_back(nuMRef.at(0));
        collectionRef.push_back(tauPRef.at(0));
        sort(piPRef.begin(),piPRef.end(), ACTools::cmpPt<reco::GenParticleRef>);
        for (unsigned int i=0; i!=piPRef.size(); i++) collectionRef.push_back(piPRef.at(i));
        collectionRef.push_back(nuPRef.at(0));
        sort(quarksRef.begin(),quarksRef.end(), ACTools::cmpPz<const reco::GenParticleRef>);
        for (unsigned int i=0; i!=quarksRef.size(); i++)    collectionRef.push_back(quarksRef.at(i));

        // printf("nuM = (%8.6f, %8.6f, %8.6f)\n", nuMRef.at(0)->vx(), nuMRef.at(0)->vy(), nuMRef.at(0)->vz());
        // printf("nuP = (%8.6f, %8.6f, %8.6f)\n", nuPRef.at(0)->vx(), nuPRef.at(0)->vy(), nuPRef.at(0)->vz());



    } else {
        std::cout<<"evt"<<iEvent.id().event()<<" GenSelector::storeEventVBFH: H->tau+tau->2x3prong not found\t";
        // std::cout<<"#of tau "<<tauCand->size()<<"#of nu "<<nuCand->size()<<"#of pi "<<piCand->size()<<std::endl;
        // std::cout<<"evtType"<<std::endl;
        std::cout<<"sizes: nus "<<nuM.size()<<nuP.size()<<", pis "<<piM.size()<<piP.size()<<std::endl;
        std::cout<<"gen-";
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
