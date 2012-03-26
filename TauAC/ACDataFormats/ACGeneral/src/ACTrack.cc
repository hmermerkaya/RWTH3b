#include "../interface/ACTrack.h"

ClassImp(ACTrack)

ACTrack::ACTrack():
ACCandidate() {
    chi2_ = -1.;
    ndof_ = -1.;
    vertex_.SetXYZ(-1000,-1000,-1000);
}
ACTrack::ACTrack(const TLorentzVector & p4, float charge, const std::string & name, float chi2, float ndof, const TVector3 & vertex):
ACCandidate(p4, charge, name) {
    chi2_ = chi2;
    ndof_ = ndof;
    vertex_ = vertex;
}

const TVector3 & ACTrack::vertex() const { return vertex_; }

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
double ACTrack::vx() const { return vertex().X(); }
double ACTrack::vy() const { return vertex().Y(); }
double ACTrack::vz() const { return vertex().Z(); }
