// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      KinematicTauProducer
// 
/**
 This framework module copies an existing PFTauCollection and modifies the tau's parameters according to a kinematical refit of its decay.
 A particular decay mode has to be assumed. New discriminators regarding the fits quality are provided.
 Part of the KinematicTau package.

 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#ifndef KinematicTauProducer_h
#define KinematicTauProducer_h


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"




class KinematicTauProducer : public edm::EDProducer {
public:
  explicit KinematicTauProducer(const edm::ParameterSet&);
  ~KinematicTauProducer();
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  //Execute the kinematic fit and in case of success modify the taus parameters
  bool select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,reco::RecoChargedCandidateCollection & daughterCollection,const edm::EventSetup& iSetup);
  //combine a discriminator of important quality cuts of a refitted tau decay
  bool dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr, const int & fitStatus, SelectedKinematicDecay &KFTau);
  int  saveKinParticles(const KinematicTauCreator * kinTauCrtr, SelectedKinematicDecayCollection &refitDecays, std::map<std::string, bool> tauDiscriminators);
  void saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection);
  void correctReferences(SelectedKinematicDecayCollection & selected, const edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands);
  void setMissingQualityCriteria(SelectedKinematicDecayCollection &decay, const KinematicTauCreator * kinTauCrtr);

  
  const edm::ParameterSet fitParameters_;
  edm::Event * iEvent_;
  edm::InputTag KinematicTauCandTag_;
  unsigned int cnt_, cntFound_;
	
};
#endif
