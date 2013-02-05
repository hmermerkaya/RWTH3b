#include "../interface/ACEventGlobals.h"

ClassImp(ACEventGlobals)

ACEventGlobals::ACEventGlobals(const ACMET & pfMET, const ACMET & pfType1CorrectedMET, const ACMET & tcMET) {
    pfMET_ = pfMET;
    pfType1CorrectedMET_ = pfType1CorrectedMET;
    tcMET_ = tcMET;
}
