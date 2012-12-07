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

#include "TTree.h"
#include "TFile.h"



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
  bool FitKinematicTauCandidate(SelectedKinematicDecay &KFTau,std::vector<reco::TrackRef> &usedTracks, edm::ESHandle<TransientTrackBuilder> &transTrackBuilder_,edm::Handle<reco::GenParticleCollection> &genParticles);
  //combine a discriminator of important quality cuts of a refitted tau decay
  bool dicriminatorByKinematicFitQuality(unsigned int &ambiguity,const KinematicTauCreator *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau);
  int  saveKinParticles(unsigned int &ambiguity,const KinematicTauCreator * kinTauCreator, SelectedKinematicDecay &KFTau);
  void saveSelectedTracks(const std::vector<reco::TrackRef> usedTracks, reco::RecoChargedCandidateCollection & daughterCollection);
  void correctReferences(SelectedKinematicDecayCollection & selected, const edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands);
  void fillTree(std::vector<double> &QCVar);
  double VertexRotationAndSignificance(TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,
                                       TVector3 &tauFlghtDirNoCorr,
                                       reco::Vertex &pVtx, TLorentzVector &lorentzA1,
                                       TVector3 &tauFlghtDir,double &theta0, double &thetaMax);
  
  const edm::ParameterSet fitParameters_;
  edm::Event * iEvent_;
  edm::InputTag primVtxTag_,KinematicTauCandTag_;
  std::vector<std::string> VertexTags_;
  std::vector<std::string> TauVtxList_;
  unsigned int cnt_, cntFound_;
  edm::InputTag gensrc_;
  unsigned int minTau_, minVtxTracks_;
  double etacut_, sigcut_;

  // BDT variables
  TFile *output;
  TTree *output_tree;

  double BDT_vtxSignPVRotSV;
  double BDT_vtxSignPVRotPVRed;
  double BDT_a1Mass;
  double BDT_energyTFraction;
  double BDT_chiSquared;

};
#endif
