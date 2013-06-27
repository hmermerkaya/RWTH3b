// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      ThreeProngInputSelector_Step2
// 
/**
 * This framework module creates the appropriate input for the KinematicTauProducer.
 * In this case the tau decay tau->3pi+nu is assumed.
 * It creates a collection of reco::TrackRefVector of size 3 for each tau candidate and refits the primary vertex without tracks assigned already to the tau decay.
 * Part of the KinematicTau package.
 *
 * @author Lars Perchalla, Philip Sauerland
 * @date 2012
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Feb  18 13:25:12 CEST 2010
//
//


#include "RecoTauTag/KinematicTau/interface/KinematicTauTools.h"

// system include files
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

#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

#include <TLorentzVector.h>
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class ThreeProngInputSelector_Step2 : public edm::EDFilter, protected KinematicTauTools {
public:
  ThreeProngInputSelector_Step2(const edm::ParameterSet&);
  virtual ~ThreeProngInputSelector_Step2();
	
private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool select(std::vector<reco::TrackRefVector> & selected, reco::PFTauRefVector & taurefs, reco::VertexCollection & primaryVertex);
  
  edm::Event * iEvent_;
  edm::ParameterSet iConfig_;
  edm::InputTag threeProngCollectionTag_, selectedTauCandidatesTag_, primVtxTag_;
  unsigned int cnt_, cntFound_, minTau_, minVtxTracks_;
  double maxChi2ndf_;    

};
