#include "../interface/ACEventInfo.h"

ClassImp(ACEventInfo)

ACEventInfo::ACEventInfo(int eventID, int lumiID, int runID, std::string type) {
    eventID_ = eventID;
    lumiID_ = lumiID;
    runID_ = runID;
    type_ = type;
}
ACEventInfo::ACEventInfo() {
    ACEventInfo(-1, -1, -1, "");
}
