#ifndef ACMuonConverter_h
#define ACMuonConverter_h
/**
 Particle helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACMuon.h"
#include "ACFrameworkModules/Common/interface/ConversionLog.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

class ACMuonConverter : public ACMuon {
public:
    virtual std::string classname() const { return "ACMuonConverter"; }
    /// constructor from reco::Muon to initiate a ACMuon
    ACMuonConverter(const reco::Muon& particle, bool isRecommendedMuon = false, const TrackConversionLog * conversionLogTrack = 0, ACGenParticleRef genRef = ACGenParticleRef());
    virtual ~ACMuonConverter() {};
};

#endif
