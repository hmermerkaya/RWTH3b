#include "../interface/ACJet.h"

ClassImp(ACJet)

ACJet::ACJet() {
    p4_.SetXYZM(-1000,-1000,-1000,-1000);
    charge_ = -1000.0;
    multiplicity_ = -1;
}
ACJet::ACJet(const TLorentzVector& p4, int multiplicity, int charge) {
    p4_ = p4;
    multiplicity_ = multiplicity;
    charge_ = charge;
}
