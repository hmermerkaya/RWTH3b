#include "../interface/ACTriggerObjectConverter.h"

ACTriggerObjectConverter::ACTriggerObjectConverter(const trigger::TriggerObject & triggerObject) {
    pdgID_ = triggerObject.id();
    pt_ = triggerObject.pt();
    eta_ = triggerObject.eta();
    phi_ = triggerObject.phi();
    mass_ = triggerObject.mass();
}
