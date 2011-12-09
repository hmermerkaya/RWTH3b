#include "ACFrameworkModules/Common/interface/PFTauMatching.h"

PFTauMatching::PFTauMatching() {
}

void PFTauMatching::nextEvt() {
    conversionLog_.clear();
}
void PFTauMatching::logConversion(const reco::PFTauRef& onlineRef, const ACPFTauRef & offlineRef) {
    conversionLog_.insert(std::make_pair(onlineRef, offlineRef));
}
const ACPFTauRef * PFTauMatching::getConverted(const reco::PFTauRef& onlineRef) const {
    for (std::map<reco::PFTauRef,ACPFTauRef>::const_iterator entry = conversionLog_.begin();  entry != conversionLog_.end(); ++entry) {
        if (onlineRef == entry->first) return &entry->second;
    }
    return 0;
}
