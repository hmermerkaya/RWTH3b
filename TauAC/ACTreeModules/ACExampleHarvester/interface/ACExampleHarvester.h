#ifndef ACExampleHarvester_h
#define ACExampleHarvester_h
/**
 Example of a fully working harvester derived from the appropiate base class.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TStopwatch.h>

#include "../../ACHarvester/interface/ACHarvester.h"

class ACExampleHarvester : public ACHarvester {
public:
    ACExampleHarvester(const std::string & outpath = ".", bool modifyInputFiles = false);

protected:
//    /// compare histograms of several dataset one-by-one (can be replaced by user-defined function)
//    void compareHistograms(TH1 * hist1, TH1 * hist2, const ACDataset * dataset1, const ACDataset * dataset2, const std::string & storagePath, const std::string & path);
//    /// persistantly modify dataset, only called if modifyInputFiles set to true. returns whether or not changes need to be stored
//    bool modifyDataset(ACDataset * dataset);
};

#endif
