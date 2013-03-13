// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      ThreeProngInputSelector_Step1
// 
/**
 The ThreeProngInputSelector_Step1 selects tracks from within a PFTaus signal cone and stores combinations of these according to the required number of tracks defined by the user.
 
 This framework module returns a boolean whether the minimum number of taus was found with enough daughters.
 For convenience a reference to the taus is also stored in case of success.
 Part of the KinematicTau package.
 
 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#ifndef ThreeProngInputSelector_Step1_h
#define ThreeProngInputSelector_Step1_h


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

#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TString.h"

class ThreeProngInputSelector_Step1 : public edm::EDProducer {
public:
  explicit ThreeProngInputSelector_Step1(const edm::ParameterSet&);
  ~ThreeProngInputSelector_Step1();
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool select(std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates,std::vector<reco::TrackCollection> &NonTauTracksLists_);
  reco::TrackRefVector getPFTauDaughters(reco::PFTauRef &PFTau);
  bool GetNonTauTracks(edm::Event *iEvent,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks, std::vector<reco::TrackRef> &tautracks);
  bool GetNonTauTracksFromVertex(SelectedKinematicDecay cand,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks);
  bool sumCharge(const std::vector<reco::TrackRef> & input);
  std::vector<std::vector<reco::TrackRef> > choose3Prongs(std::vector<reco::TrackRef> & input);
  std::vector<std::vector<reco::TrackRef> > permuteCombinations(const std::vector<reco::TrackRef> & vect);
  double getInvariantMass(std::vector<reco::TrackRef> &tracks);
  template <typename T> static bool cmpPt(const T & a, const T & b);
  
  edm::Event * iEvent_;
  std::string tauType_;
  edm::InputTag trkCollectionTag_,vtxtrackCollectionTag_, primVtx_;
  edm::ProductID trkCollectionID_;
  unsigned int minTracks_, minTau_, nTauPerVtx_, cnt_, cntFound_;
  double minTauPt_,TauEtaCut_;
  std::vector<std::string> TauVtxList_;
  ParticleMassHelper PMH;

};

#endif
