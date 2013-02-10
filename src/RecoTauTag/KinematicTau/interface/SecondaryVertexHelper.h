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
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class SecondaryVertexHelper {
public:
  SecondaryVertexHelper(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder,const SelectedKinematicDecay &KTau);
  ~SecondaryVertexHelper();
  
  bool                              hasSecondaryVertex(){return hasSecondaryVertex_;}
  TransientVertex                   InitialSecondaryVertex(){return tmpVtx_;}
  std::vector<reco::TransientTrack> InitialRefittedTracks(){return trks_;}
  TLorentzVector                    Initial_a1_p4(){return a1_p4_;}
  std::vector<TLorentzVector>       Initial_pions(){return pions_;}

private:
  bool checkSecVtx(std::vector<reco::TransientTrack> trkVct, TransientVertex & transVtx);

  bool hasSecondaryVertex_;
  TransientVertex tmpVtx_;
  std::vector<reco::TransientTrack> trks_;
  TLorentzVector a1_p4_;
  std::vector<TLorentzVector> pions_;

};

#endif
