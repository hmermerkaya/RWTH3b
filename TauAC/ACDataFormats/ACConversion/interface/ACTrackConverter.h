#ifndef ACTrackConverter_h
#define ACTrackConverter_h
/**
 Track helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACTrack.h"

#include <DataFormats/TrackReco/interface/Track.h>

class ACTrackConverter : public ACTrack {
public:
    virtual std::string classname() const { return "ACTrackConverter"; }
    /// constructor from reco::Track to initiate a ACTrack
    ACTrackConverter(const reco::Track& track);
    virtual ~ACTrackConverter() {};
};

#endif
