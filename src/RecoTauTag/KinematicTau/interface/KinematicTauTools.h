// -*- C++ -*-
//
// Package: RecoTauTag/KinematicTau    
// Class: KinematicTauTools     
//
// Original Author:  Ian Nugent
//         Created:  Mon July 17 2012
//                                                                                                                                                                                                                                           
#ifndef KinematicTauTools_h
#define KinematicTauTools_h
// system include files                                                                                                                                                                                              
#include <memory>
#include <fstream>

// user include files                                                                                                                                                                                                
#include <memory>
#include <fstream>

// user include files                                                                                                                                                                                           
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TLorentzVector.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

typedef std::vector<std::vector<std::vector<reco::TrackRef> > > vVVTrackRef;
typedef std::vector<std::vector<reco::TrackRef> > vVTrackRef;

class KinematicTauTools {

public:
  // constructor and Destructor
  KinematicTauTools();
  virtual ~KinematicTauTools();

  bool sumCharge(const std::vector<reco::TrackRef> & input);
  vVTrackRef choose3Prongs(std::vector<reco::TrackRef> & input);
  double VertexRotationAndSignificance(const std::vector<reco::TrackRef> & input,TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,reco::Vertex &pVtx,TLorentzVector &lorentzA1, TVector3 &tauFlghtDir, double &theta0, double &thetaMax);
  bool choose3bestTracks(std::vector<reco::TrackRef> & input, reco::Vertex & pVtx);
  bool choose3bestTracks(std::vector<reco::TrackRefVector> & selected, std::vector<std::vector<reco::TrackRef> > combis, const reco::Vertex & pVtx);
  bool removeDuplicateTriplets(const std::vector<reco::TrackRef> & duplicateTracks, vVVTrackRef & threeProngCombis, vVVTrackRef::iterator & candidates, vVTrackRef::iterator & triplets);
  bool checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx);
  std::vector<reco::TransientTrack> convToTransTrck(std::vector<reco::TrackRef> &input);

  //template classes
  template <typename T> std::vector<std::vector<T> > permuteCombinations(const std::vector<T> & vect);
  template <typename T> double getInvariantMass(const T & tracks);
  template <class T> TLorentzVector getSumTLorentzVec(const T& tracks, const double massConstraint);
  template <typename T> static bool cmpPt(const T & a, const T & b);
  template <typename S, typename T> static bool pairSecond(const std::pair<S, T> &a, const std::pair<S, T> &b);
  template <typename T> static bool cmpChi2(const T &a, const T &b);

  inline double Get_piMass(){return piMass;}
  inline double Get_tauMass(){return tauMass;}

  void Set_TransientTrackBuilder(edm::ESHandle<TransientTrackBuilder>  transTrackBuilder){transientTrackBuilder_=transTrackBuilder;}
  bool GetNonTauTracks(edm::Event *iEvent,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks, std::vector<reco::TrackRef> &tautracks);

protected:
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder_;

private:
  static const double piMass;
  static const double tauMass;
  
};

#endif
