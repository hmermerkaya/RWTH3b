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


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

#include <TLorentzVector.h>
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class ThreeProngInputSelector_Step2 : public edm::EDProducer {
public:
  ThreeProngInputSelector_Step2(const edm::ParameterSet&);
  virtual ~ThreeProngInputSelector_Step2();
	
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool select(std::vector<SelectedKinematicDecay> &PreKinematicDecaysStep2_,const edm::EventSetup& iSetup);
  double VertexRotationAndSignificance(TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,
				       TVector3 &tauFlghtDirNoCorr,
				       reco::Vertex &pVtx, TLorentzVector &lorentzA1,
				       TVector3 &tauFlghtDir,double &theta0, double &thetaMax);

  
  edm::Event * iEvent_;
  edm::ParameterSet iConfig_;
  edm::InputTag primVtxTag_,KinematicTauCandTag_;
  unsigned int cnt_, cntFound_, minTau_, minVtxTracks_;
  double maxChi2ndf_;    
};
