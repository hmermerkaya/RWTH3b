// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      KinematicTauCreator
// 
/**
 This is a pure abstract class providing the interface to the KinematicFit. Use derived classes e.g. ThreeProngTauCreator to implement the create() function.
 Part of the KinematicTau package.

 @author Lars Perchalla, Philip Sauerland
 @date 2010
 */

#ifndef KinematicTauCreator_h
#define KinematicTauCreator_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

class KinematicTauCreator {
public:
  KinematicTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder);
  KinematicTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const edm::ParameterSet& cfg);
  virtual ~KinematicTauCreator();
  
  virtual int create(unsigned int& ambiguity,SelectedKinematicDecay &KFTau) = 0;

  reco::PFTau getPFTau() const;
  reco::PFTau getKinematicTau() const;
  std::vector<math::XYZTLorentzVector> getRefittedChargedDaughters() const;
  std::vector<math::XYZTLorentzVector> getRefittedNeutralDaughters() const;
  std::vector<reco::TrackRef> getSelectedTracks() const;
  RefCountedKinematicTree getKinematicTree() const;
  KinematicConstrainedVertexFitter * getFitter() const {return kcvFitter_;}
  //KinematicParticleVertexFitter * getFitter() const {return kcvFitter_;}
  reco::Vertex getModifiedPrimaryVertex() const;
  
  float chi2() const ;
  virtual int ndf() const = 0;
  
protected:
  KinematicConstrainedVertexFitter *kcvFitter_;
  //KinematicParticleVertexFitter *kcvFitter_;
  RefCountedKinematicTree kinTree_;
  std::vector<reco::TrackRef> selectedTracks_;
  reco::Vertex modifiedPV_;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder_;
};

#endif
