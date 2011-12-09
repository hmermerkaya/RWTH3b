/**
 Class to log the association between PFTauRefs and ACPFTauRefs.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <map>

#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "ACDataFormats/ACGeneral/interface/ACRef.h"

class PFTauMatching {
public:
    PFTauMatching();
    virtual ~PFTauMatching() {};

    /// trigger new event. cleans up counters and references.
    void nextEvt();
    /// store the association between an offline ref and the corresponding framework ref from which it was created from.
    void logConversion(const reco::PFTauRef& onlineRef, const ACPFTauRef & offlineRef);
    /// retrieve a reference to the converted offline object from a framework reference
    const ACPFTauRef * getConverted(const reco::PFTauRef & onlineRef) const;

private:
    /// association between an offline ref and the corresponding framework ref from which it was created from.
    std::map<reco::PFTauRef,ACPFTauRef> conversionLog_;
};
