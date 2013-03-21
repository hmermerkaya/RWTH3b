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
#include "RecoTauTag/KinematicTau/interface/FitSequencer.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DecisionTree.h"


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
  bool select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,const edm::EventSetup& iSetup);
  bool FitKinematicTauCandidate(SelectedKinematicDecay &KFTau, edm::ESHandle<TransientTrackBuilder> &transTrackBuilder_,edm::Handle<reco::GenParticleCollection> &genParticles);
  //combine a discriminator of important quality cuts of a refitted tau decay
  bool dicriminatorByKinematicFitQuality(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau);
  int  saveKinParticles(unsigned int &ambiguity,FitSequencer * kinTauCreator, SelectedKinematicDecay &KFTau);
  double VertexRotationAndSignificance(TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,
                                       TVector3 &tauFlghtDirNoCorr,
                                       reco::Vertex &pVtx, TLorentzVector &lorentzA1,
                                       TVector3 &tauFlghtDir,double &theta0, double &thetaMax);

  bool GetNonTauTracksFromVertex(SelectedKinematicDecay cand,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // BDT functions
  void fillTree(std::vector<double> &QCVar);
  void FillTreeForTraining(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau);
  double ReturnBDTOutput(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau);
  
  const edm::ParameterSet fitParameters_;
  edm::Event * iEvent_;
  edm::InputTag primVtxTag_,KinematicTauCandTag_;
  edm::InputTag trkCollectionTag_;
  edm::InputTag beamSpot_;
  unsigned int cnt_, cntFound_;
  std::vector<unsigned int> cntSVFound_,cntSVQC_,cntLCFit_,cntLCFitQC_;
  edm::InputTag gensrc_;
  unsigned int minTau_, minVtxTracks_;
  double etacut_, sigcut_;
  bool useTrackHelixFit_;
  bool do_BDTTrain_;
  bool do_BDTComp_;

  TMVA::Reader *reader;

  std::string BDTweightFileMinus_;
  std::string BDTweightFilePlus_;
  std::string BDTweightFileZero_;

  float fracMins;
  float a1MassMins;
  float ProbMins3;
  float iterMins;
  float PVSVMins;

  float fracPlus;
  float a1MassPlus;
  float ProbPlus3;
  float iterPlus;
  float PVSVPlus;

  float fracZero;
  float a1MassZero;
  float ProbZero3;
  float iterZero;
  float PVSVZero;

  // BDT variables
  TFile *output;
  TTree *output_tree;

  std::vector<double> BDT_vtxSignPVRotSV;
  std::vector<double> BDT_vtxSignPVRotPVRed;
  std::vector<double> BDT_a1Mass;
  std::vector<double> BDT_energyTFraction;
  std::vector<double> BDT_chiSquared;
  std::vector<double> BDT_iterations;      

};
#endif
