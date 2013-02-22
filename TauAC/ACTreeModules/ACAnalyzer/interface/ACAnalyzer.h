#ifndef ACAnalyzer_h
#define ACAnalyzer_h
/**
 Base class to communicate between the ACTreeReader module and the specific analysis.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACEvent/interface/ACEvent.h"
#include "../../ACProvenance/interface/ACProvenance.h"

class ACAnalyzer {
public:
    ACAnalyzer();
    virtual void analyze(const ACEvent & event) = 0;
    virtual void beginJob(const ACProvenance & provenance) {};
    virtual void endJob() {};
    virtual void beginRun() {};
    virtual void endRun() {};
    virtual void beginLuminosityBlock() {};
    virtual void endLuminosityBlock() {};

protected:
};

#endif
