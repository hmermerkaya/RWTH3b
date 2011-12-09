#ifndef ACTrigger_h
#define ACTrigger_h
/**
 Class containing trigger results.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCommon.h"
#include <map>

class ACTrigger : public ACCommon {
public:
    virtual std::string classname() const { return "ACTrigger"; }
    ACTrigger() {};
    ACTrigger(std::string tableName, std::map<std::string,int> triggerPathMap);
    virtual ~ACTrigger() {};

    const std::string & tableName() const;
    const std::map<std::string,int> & hltPathMap() const;

    ClassDef(ACTrigger,1);

protected:
    /// name of the trigger menu
    std::string tableName_;
    /// map of all HLT paths and their execution bit
    std::map<std::string,int> hltPathMap_;
};

typedef ACTrigger* pACTrigger;

#endif
