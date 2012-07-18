#ifndef ConversionLog_h
#define ConversionLog_h
/**
 Class to log the association between framework references and ACRefs.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <map>

#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "ACDataFormats/ACGeneral/interface/ACRef.h"

template <class S, class T> class ConversionLog {
public:
    ConversionLog() {};
    virtual ~ConversionLog() {};

    /// trigger new event. cleans up counters and references.
    void nextEvt();
    /// store the association between an offline ref and the corresponding framework ref from which it was created from.
    void logConversion(const S & onlineRef, const T & offlineRef);
    /// retrieve a reference to the converted offline object from a framework reference
    const T * getConverted(const S & onlineRef) const;

private:
    /// association between an offline ref and the corresponding framework ref from which it was created from.
    std::map<S, T> conversionLog_;
};

template<class S, class T> void ConversionLog<S, T>::nextEvt() {
    conversionLog_.clear();
}
template<class S, class T> void ConversionLog<S, T>::logConversion(const S & onlineRef, const T & offlineRef) {
    conversionLog_.insert(std::make_pair(onlineRef, offlineRef));
}
template<class S, class T> const T * ConversionLog<S, T>::getConverted(const S & onlineRef) const {
    for (typename std::map<S, T>::const_iterator entry = conversionLog_.begin();  entry != conversionLog_.end(); ++entry) {
        if (onlineRef == entry->first) return &entry->second;
    }
    return 0;
}

typedef ConversionLog<reco::PFTauRef, ACPFTauRef> PFTauConversionLog;
typedef ConversionLog<reco::PFJetRef, ACJetRef> PFJetConversionLog;
typedef ConversionLog<reco::TrackRef, ACTrackRef> TrackConversionLog;

#endif
