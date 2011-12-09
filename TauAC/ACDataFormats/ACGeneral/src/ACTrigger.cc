#include "../interface/ACTrigger.h"

ClassImp(ACTrigger)

ACTrigger::ACTrigger(std::string tableName, std::map<std::string,int> hltPathMap) {
    tableName_ = tableName;
    hltPathMap_ = hltPathMap;
}
const std::string & ACTrigger::tableName() const {
    return tableName_;
}
const std::map<std::string,int> & ACTrigger::hltPathMap() const {
    return hltPathMap_;
}
