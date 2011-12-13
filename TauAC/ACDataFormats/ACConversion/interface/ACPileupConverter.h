#ifndef ACPileupConverter_h
#define ACPileupConverter_h
/**
 Pileup helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACPileupInfo.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class ACPileupConverter : public ACPileupInfo {
public:
    virtual std::string classname() const { return "ACPileupConverter"; }
    /// constructor from PileupSummaryInfo to initiate an ACPileupInfo
    ACPileupConverter(const PileupSummaryInfo & pileup);
    virtual ~ACPileupConverter() {};
};

#endif
