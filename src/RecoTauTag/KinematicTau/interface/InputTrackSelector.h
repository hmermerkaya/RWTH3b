// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      InputTrackSelector
// 
/**
 The InputTrackSelector selects tracks from within a PFTaus signal cone and stores combinations of these according to the required number of tracks defined by the user.
 
 This framework module returns a boolean whether the minimum number of taus was found with enough daughters.
 For convenience a reference to the taus is also stored in case of success.
 Part of the KinematicTau package.
 
 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#ifndef InputTrackSelector_h
#define InputTrackSelector_h


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

#include  "RecoTauTag/KinematicTau/interface/KinematicTauTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TString.h"

class InputTrackSelector : public edm::EDProducer, protected KinematicTauTools {
public:
  explicit InputTrackSelector(const edm::ParameterSet&);
  ~InputTrackSelector();
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool select(std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates,std::vector<reco::TrackCollection> &NonTauTracksLists_);
  reco::TrackRefVector getPFTauDaughters(reco::PFTauRef &PFTau);
  
  edm::Event * iEvent_;
  std::string tauType_;
  edm::InputTag trkCollectionTag_,vtxtrackCollectionTag_, primVtx_;
  edm::ProductID trkCollectionID_;
  unsigned int minTracks_, minTau_, nTauPerVtx_, cnt_, cntFound_;
  double minTauPt_;
  std::vector<std::string> TauVtxList_;
};

#endif
