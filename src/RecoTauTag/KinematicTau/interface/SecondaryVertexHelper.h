// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      SecondaryVertexHelper
// 
/**
   Helper Class to get secondary vertex information
 
 @author Ian Nugent
 @date 2012
 */

#ifndef SecondaryVertexHelper_h
#define SecondaryVertexHelper_h

// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include  "RecoTauTag/KinematicTau/interface/KinematicTauTools.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class SecondaryVertexHelper : protected KinematicTauTools {
public:
  SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder,const SelectedKinematicDecay &KTau);
  ~SecondaryVertexHelper();
  
  bool hasSecondaryVertex(){return hasSecondaryVertex_;}
  TransientVertex SecondaryVertex(){return tmpVtx_;}
  std::vector<reco::TransientTrack> RefittedTracks(){return trks_;}

private:
  bool hasSecondaryVertex_;
  TransientVertex tmpVtx_;
  std::vector<reco::TransientTrack> trks_;
};

#endif
