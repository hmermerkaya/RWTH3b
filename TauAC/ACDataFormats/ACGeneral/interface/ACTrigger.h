#ifndef ACTrigger_h
#define ACTrigger_h
/**
 Class containing trigger results.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <map>

#include "ACCommon.h"
#include "ACRef.h"

class ACTrigger : public ACCommon {
public:
    virtual std::string classname() const { return "ACTrigger"; }
    ACTrigger() {};
    /// basic c'tor, no trigger object addons
    ACTrigger(const std::string & tableName, const std::map<std::string, int> & hltPathMap);
    /// advanced c'tor storing also trigger object addons
    ACTrigger(const std::string & tableName, const std::map<std::string, int> & hltPathMap, const std::map<std::string, std::vector<std::string> > & moduleNamesPerPath, const std::map<std::string, std::vector<ACRef<ACTriggerObject> > > & objectsPerModule);
    virtual ~ACTrigger() {};

    const std::string & tableName() const { return tableName_; }
    const std::map<std::string,int> & hltPathMap() const { return hltPathMap_; }
    const std::map<std::string, std::vector<std::string> > & moduleNamesPerPath() const { return moduleNamesPerPath_; }
    const std::map<std::string, std::vector<ACRef<ACTriggerObject> > > & objectsPerModule() const { return objectsPerModule_; }

    ClassDef(ACTrigger,1);

protected:
    /// name of the trigger menu
    std::string tableName_;
    /// map of all HLT paths and their execution bit
    std::map<std::string, int> hltPathMap_;
    /// map of all HLT paths and their module names
    std::map<std::string, std::vector<std::string> > moduleNamesPerPath_;
    /// map of all module names and their objects (as ref pointing to branch storing a list of all unique trigger objects)
    std::map<std::string, std::vector<ACRef<ACTriggerObject> > > objectsPerModule_;
};

typedef ACTrigger * pACTrigger;

#endif
