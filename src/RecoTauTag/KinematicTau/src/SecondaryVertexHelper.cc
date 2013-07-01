#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"

#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>
#include <iostream>

SecondaryVertexHelper::SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const SelectedKinematicDecay &KTau):
  hasSecondaryVertex_(false),
  a1_p4_(0,0,0,0)
{
  ParticleMassHelper PMH;
  std::vector<reco::TrackRef> input=KTau.InitialTrackTriplet();
  for(unsigned int i=0; i<input.size();i++){
    pions_.push_back(TLorentzVector(input.at(i)->px(),input.at(i)->py(),input.at(i)->pz(),sqrt(pow(input.at(i)->p(),2.0)+pow(PMH.Get_piMass(),2.0))));
  }
  // build transient tracks
  for (std::vector<reco::TrackRef>::const_iterator iter=input.begin(); iter!=input.end(); ++iter) {
    trks_.push_back( transTrackBuilder->build( *iter ) );
  }
  // get vertex
  if (checkSecVtx(trks_,tmpVtx_)){
    hasSecondaryVertex_=true;
    // compute a1 vector at secondary vertex
    for(unsigned int i=0; i<trks_.size();i++){
      TrajectoryStateClosestToPoint TPCTP=trks_.at(i).trajectoryStateClosestToPoint(tmpVtx_.position());
      TVector3 pi_mom;
      pi_mom.SetPtThetaPhi(TPCTP.pt(),TPCTP.perigeeParameters().theta(),TPCTP.perigeeParameters().phi());
      pions_.at(i).SetVectM(pi_mom,PMH.Get_piMass());
    }
  }
  a1_p4_.SetPxPyPzE(0.0,0.0,0.0,0.0);
  for(unsigned int i=0; i<input.size();i++){
    a1_p4_+=pions_.at(i);
  }
}


SecondaryVertexHelper::~SecondaryVertexHelper(){

}



bool SecondaryVertexHelper::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx){
  if(trkVct.size()<2){
    LogTrace("ThreeProngTauCreator")<<"Can't check SecVertex: Only "<<trkVct.size()<<" Tracks.";
    return false;
  }else{
    bool useAdaptive = false;
    if(useAdaptive){
      AdaptiveVertexFitter avf;
      avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception                                                                                                                                            
      try{
        transVtx = avf.vertex(trkVct); //AdaptiveVertexFitter                                                                                                                                                                                
      }catch(...){
        //LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
        return false;
      }
    }else{
      KalmanVertexFitter kvf(true);
      try{
        transVtx = kvf.vertex(trkVct); //KalmanVertexFitter                                                                                                                                                                                  
      }catch(...){
        //LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
        return false;
      }
    }
    if(!transVtx.isValid()) LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex not valid.";
    if(!useAdaptive){
      if(!transVtx.hasRefittedTracks()){
        //LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has 0 refitted tracks.";
        return false;
      }else if(transVtx.refittedTracks().size()!=trkVct.size()){
        //LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has only "<<transVtx.refittedTracks().size()<<" refitted of "<<trkVct.size()<<" initial tracks.";
        return false;
      }
    }

    return transVtx.isValid();
  }
}

