#ifndef ACExampleAnalyzer_h
#define ACExampleAnalyzer_h
/**
 Example of a fully working analyzer derived from the appropiate base class.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACAnalyzer/interface/ACAnalyzer.h"
#include "../../ACTreeReader/interface/ACTreeReader.h"

class ACExampleAnalyzer : public ACAnalyzer {
public:
    ACExampleAnalyzer();
    virtual void analyze(const ACEvent & event);
    virtual void beginRun();
    virtual void endRun();
    virtual void beginLuminosityBlock();
    virtual void endLuminosityBlock();
    virtual void beginJob(const ACProvenance & provenance);
    virtual void endJob();

protected:
};

#endif
