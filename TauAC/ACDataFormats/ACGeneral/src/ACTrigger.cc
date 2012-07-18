#include "../interface/ACTrigger.h"

ClassImp(ACTrigger)

ACTrigger::ACTrigger(const std::string & tableName, const std::map<std::string, int> & hltPathMap) {
    tableName_ = tableName;
    hltPathMap_ = hltPathMap;
    moduleNamesPerPath_.clear();
    objectsPerModule_.clear();
}
ACTrigger::ACTrigger(const std::string & tableName, const std::map<std::string, int> & hltPathMap, const std::map<std::string, std::vector<std::string> > & moduleNamesPerPath, const std::map<std::string, std::vector<ACRef<ACTriggerObject> > > & objectsPerModule){
    tableName_ = tableName;
    hltPathMap_ = hltPathMap;
    moduleNamesPerPath_ = moduleNamesPerPath;
    objectsPerModule_ = objectsPerModule;
}
