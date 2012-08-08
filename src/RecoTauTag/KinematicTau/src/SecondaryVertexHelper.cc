#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"


SecondaryVertexHelper::SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const SelectedKinematicDecay &KTau):
  KinematicTauTools(),
  hasSecondaryVertex_(false)
{
  Set_TransientTrackBuilder(transTrackBuilder);
  std::vector<reco::TrackRef> input=KTau.InitalTrackTriplet();
  for(unsigned int i=0; i<input.size();i++){
    pions_.push_back(TLorentzVector(input.at(i)->px(),input.at(i)->py(),input.at(i)->pz(),sqrt(pow(input.at(i)->p(),2.0)+pow(Get_piMass(),2.0))));
  }
  //double massA1 = getInvariantMass(input);
  //a1_p4_ = getSumTLorentzVec(input, massA1);
  trks_ = convToTransTrck(input);
  if (checkSecVtx(trks_,tmpVtx_))hasSecondaryVertex_=true;
}

SecondaryVertexHelper::~SecondaryVertexHelper(){

}
