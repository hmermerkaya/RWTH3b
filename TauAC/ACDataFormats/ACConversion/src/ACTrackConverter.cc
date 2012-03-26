#include "../interface/ACTrackConverter.h"

ACTrackConverter::ACTrackConverter(const reco::Track& track) {
    p4_ = TLorentzVector(track.px(), track.py(), track.pz(), 0.);
    charge_ = track.charge();
    name_ = "track";

    chi2_ = track.chi2();
    ndof_ = track.ndof();
    vertex_ = TVector3(track.vertex().X(), track.vertex().Y(), track.vertex().Z());
}
