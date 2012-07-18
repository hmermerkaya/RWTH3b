#ifndef ACBeamSpotConverter_h
#define ACBeamSpotConverter_h
/**
 Helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla
 @date 2012
 */

#include "../../ACGeneral/interface/ACBeamSpot.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class ACBeamSpotConverter : public ACBeamSpot {
public:
    virtual std::string classname() const { return "ACBeamSpotConverter"; }
    /// constructor from reco::Vertex to initiate an ACBeamSpot
    ACBeamSpotConverter(const reco::BeamSpot& beamSpot);
    virtual ~ACBeamSpotConverter() {};
};

#endif
