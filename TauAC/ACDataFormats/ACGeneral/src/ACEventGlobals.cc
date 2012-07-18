#include "../interface/ACEventGlobals.h"

ClassImp(ACEventGlobals)

ACEventGlobals::ACEventGlobals(const TVector3& pfMET, const TVector3& pfType1CorrectedMET, const TVector3& tcMET, double pfSumET, double pfType1CorrectedSumET, double tcSumET) {
    pfMET_ = pfMET;
    pfType1CorrectedMET_ = pfType1CorrectedMET;
    tcMET_ = tcMET;
    pfSumET_ = pfSumET;
    pfType1CorrectedSumET_ = pfType1CorrectedSumET;
    tcSumET_ = tcSumET;
}
ACEventGlobals::ACEventGlobals() {
    ACEventGlobals(TVector3(-1000,-1000,-1000), TVector3(-1000,-1000,-1000), TVector3(-1000,-1000,-1000), -1000.0, -1000.0, -1000.0);
}
double ACEventGlobals::pfMETSignificance() const {
    if (pfSumET_ > 0.0) {
        return pfMET_.Mag()/sqrt(pfSumET_);
    }
    return -1.0;
}
double ACEventGlobals::tcMETSignificance() const {
    if (tcSumET_ > 0.0) {
        return tcMET_.Mag()/sqrt(tcSumET_);
    }
    return -1.0;
}
double ACEventGlobals::pfMETTopology() const {
    if (pfSumET_ > 0.0) {
        return pfMET_.Mag()/pfSumET_;
    }
    return -1.0;
}
double ACEventGlobals::tcMETTopology() const {
    if (tcSumET_ > 0.0) {
        return tcMET_.Mag()/tcSumET_;
    }
    return -1.0;
}
