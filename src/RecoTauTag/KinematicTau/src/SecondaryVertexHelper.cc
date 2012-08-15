#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"


SecondaryVertexHelper::SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const SelectedKinematicDecay &KTau):
  KinematicTauTools(),
  hasSecondaryVertex_(false),
  a1_p4_(0,0,0,0)

{
  Set_TransientTrackBuilder(transTrackBuilder);
  std::vector<reco::TrackRef> input=KTau.InitalTrackTriplet();
  for(unsigned int i=0; i<input.size();i++){
    pions_.push_back(TLorentzVector(input.at(i)->px(),input.at(i)->py(),input.at(i)->pz(),sqrt(pow(input.at(i)->p(),2.0)+pow(Get_piMass(),2.0))));
  }
  trks_ = convToTransTrck(input);
  if (checkSecVtx(trks_,tmpVtx_)){
    hasSecondaryVertex_=true;
    // compute a1 vector at secondary vertex
    for(unsigned int i=0; i<trks_.size();i++){
      TrajectoryStateClosestToPoint TPCTP=trks_.at(i).trajectoryStateClosestToPoint(tmpVtx_.position());
      TVector3 pi_mom;
      pi_mom.SetPtThetaPhi(TPCTP.pt(),TPCTP.perigeeParameters().theta(),TPCTP.perigeeParameters().phi());
      pions_.at(i).SetVectM(pi_mom,Get_piMass());
    }
  }
  for(unsigned int i=0; i<input.size();i++){
    a1_p4_+=pions_.at(i);
  }
}

SecondaryVertexHelper::~SecondaryVertexHelper(){

}
