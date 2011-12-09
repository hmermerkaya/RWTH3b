#include "../interface/ACCandidate.h"

ClassImp(ACCandidate)

ACCandidate::ACCandidate() {
    p4_.SetXYZM(-1000,-1000,-1000,-1000);
}
ACCandidate::ACCandidate(const TLorentzVector& p4, float charge) {
    p4_ = p4;
    charge_ = charge;
}
float ACCandidate::charge() const { return charge_; }
const TLorentzVector & ACCandidate::p4() const { return p4_; }

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
double ACCandidate::p() const { return p4().P(); }
double ACCandidate::energy() const { return p4().E(); }
double ACCandidate::et() const { return p4().Et(); }
double ACCandidate::mass() const { return p4().M(); }
double ACCandidate::mt() const { return p4().Mt(); }
double ACCandidate::px() const { return p4().Px(); }
double ACCandidate::py() const { return p4().Py(); }
double ACCandidate::pz() const { return p4().Pz(); }
double ACCandidate::pt() const { return p4().Pt(); }
double ACCandidate::phi() const { return p4().Phi(); }
double ACCandidate::theta() const { return p4().Theta(); }
double ACCandidate::eta() const { return p4().Eta(); }
