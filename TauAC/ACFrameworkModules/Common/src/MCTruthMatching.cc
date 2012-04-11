#include "ACFrameworkModules/Common/interface/MCTruthMatching.h"

MCTruthMatching::MCTruthMatching(edm::InputTag mapTag) {
    mapTag_ = mapTag;
    passedParticles_ = 0;
    cntMatched_ = 0;
    cntDuplicateMatched_ = 0;
}

void MCTruthMatching::nextEvt() {
    conversionLog_.clear();
    assignedRefsPerEvt_.clear();
}
void MCTruthMatching::logConversion(const reco::GenParticleRef& onlineRef, const ACGenParticleRef & offlineRef) {
    conversionLog_.insert(std::make_pair(onlineRef, offlineRef));
}
const ACGenParticleRef * MCTruthMatching::getMatching(const edm::Event& evt, const SelectedKinematicParticle * kinparticle) {
    passedParticles_++;
    bool foundProduct = evt.getByLabel(mapTag_, map_);
    if (!foundProduct) {
        printf("TrackMCTrth::getMatching: ERROR! branch %s not available. Continue processing.\n", mapTag_.label().c_str());
        return 0;
    }

    if (!kinparticle->candRef().isAvailable()) {
        printf("evt %d MCTruthMatching::getMatching: ERROR! No reco::Candidate available.\n", evt.id().event());
        return 0;
    }
    const reco::GenParticleRef & ref = (*map_)[kinparticle->candRef()];
//    printf("evt %d MCTruthMatching::getMatching: map size is %lu.\n", evt.id().event(), map_->size());
    if (!ref.isAvailable()) {
//        printf("evt %d MCTruthMatching::getMatching: No matching found.\n", evt.id().event());
        return 0;
    }
    if (ref->charge() != kinparticle->candRef()->charge()) {
        //printf("evt %d MCTruthMatching::getMatching: Bad charge! reco (charge, pt): %i, %f. gen (charge, pt, pdgID) %i, %f, %d.\n", evt.id().event(), kinparticle->candRef()->charge(), kinparticle->candRef()->pt(), ref->charge(), ref->pt(), ref->pdgId());
        return 0;
    }
    if (fabs(ref->pt() - kinparticle->candRef()->pt()) > .1*kinparticle->candRef()->pt()) {
        //printf("evt %d MCTruthMatching::getMatching: Bad momentum! reco (charge, pt): %i, %f. gen (charge, pt, pdgID) %i, %f, %d.\n", evt.id().event(), kinparticle->candRef()->charge(), kinparticle->candRef()->pt(), ref->charge(), ref->pt(), ref->pdgId());
        return 0;
    }
    if (kinparticle->ambiguity() != 2) { //particles from ambiguity are already tested
        if (!testDuplicate(ref)) {
//            printf("evt %d MCTruthMatching::getMatching: WARNING! duplicate match to gen %d.\n", evt.id().event(), ref.index());
            cntDuplicateMatched_++;
        }
    }

    cntMatched_++;

    return getConverted(ref);
}
const ACGenParticleRef * MCTruthMatching::getConverted(const reco::GenParticleRef& onlineRef) const {
    for (std::map<reco::GenParticleRef, ACGenParticleRef>::const_iterator entry = conversionLog_.begin();  entry != conversionLog_.end(); ++entry) {
        if (onlineRef == entry->first) return &entry->second;
    }
    printf("MCTruthMatching::getConverted: ERROR! Provided generator particle with pdgID %i not available. There are %lu conversions of generator particles logged.\n", onlineRef->pdgId(), conversionLog_.size());
    return 0;

}
void MCTruthMatching::increment(bool matched) {
    passedParticles_++;
    if (matched) cntMatched_++;
}
void MCTruthMatching::printOutro() {
    float ratio = 0.0;
    if (passedParticles_!=0) ratio = (float)cntMatched_/passedParticles_;
    printf("--> [MCTruthMatching] Efficiency: %d/%d = %.2f%% including %i duplicate matches.\n", cntMatched_, passedParticles_, ratio*100.0, cntDuplicateMatched_);
}
bool MCTruthMatching::testDuplicate(const reco::GenParticleRef & ref) {
    bool noDuplicateFound = true;
    for (std::vector<reco::GenParticleRef>::const_iterator iter=assignedRefsPerEvt_.begin(); iter!=assignedRefsPerEvt_.end(); ++iter) {
        if (*iter == ref) {
            noDuplicateFound = false;
            break;
        }
    }
    assignedRefsPerEvt_.push_back(ref);

    return noDuplicateFound;
}