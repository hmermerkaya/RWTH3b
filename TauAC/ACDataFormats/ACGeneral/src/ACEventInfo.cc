#include "../interface/ACEventInfo.h"

ClassImp(ACEventInfo)

ACEventInfo::ACEventInfo(int eventID, int lumiID, int runID, std::string type) {
    eventID_ = eventID;
    lumiID_ = lumiID;
    runID_ = runID;
    type_ = type;
    filterResults_.clear();
}
ACEventInfo::ACEventInfo() {
    ACEventInfo(-1, -1, -1, "");
}

void ACEventInfo::storeEDMFilterResult(const std::string & name, int value) {
    filterResults_.push_back(std::make_pair(name, value));
}
