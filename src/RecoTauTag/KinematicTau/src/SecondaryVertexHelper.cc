#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"


SecondaryVertexHelper::SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const SelectedKinematicDecay &KTau):
  KinematicTauTools(),
  hasSecondaryVertex_(false)
{
  Set_TransientTrackBuilder(transTrackBuilder);
  std::vector<reco::TrackRef> input=KTau.TrackTriplet();
  trks_ = convToTransTrck(input);
  if (checkSecVtx(trks_,tmpVtx_))hasSecondaryVertex_=true;
}

SecondaryVertexHelper::~SecondaryVertexHelper(){

}
