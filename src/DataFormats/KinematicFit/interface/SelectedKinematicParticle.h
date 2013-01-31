#ifndef DataFormats_KinematicFit_SelectedKinematicParticle_h
#define DataFormats_KinematicFit_SelectedKinematicParticle_h
// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      SelectedKinematicParticle
// 
/**
 This class stores the results from kinematically refitted particles.

 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#include <memory>
#include <algorithm>
#include <TMath.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"

class SelectedKinematicParticle {
 public:
  SelectedKinematicParticle();
  SelectedKinematicParticle(LorentzVectorParticle p,const int status, const int ambiguity, const reco::RecoChargedCandidateRef & CandRef);
  
  const int status() const;
  const int pdgid() const;
  const int charge() const;
  const unsigned int ambiguity() const;
  
  const TVectorT<double> & parameters() const;
  const TVectorT<double> & input_parameters() const;
  
  const TMatrixDSym & matrix() const;
  const TMatrixDSym & input_matrix() const;
  
  const reco::RecoChargedCandidateRef & candRef() const;
  void setCandRef(const reco::RecoChargedCandidateRef & parm);
  
  const TLorentzVector p4() const;
  const TVector3 vertex() const;
  
  void setInitialState(const TLorentzVector & momentum, const reco::Vertex & primVtx);
  
 private:
  /// status code (not in use yet)
  int status_; //
  /// particle name
  int pdgid_; //
  /// particle charge
  int charge_; //
  /// DEPRECATED! ambiguity counter
  int ambiguity_; //
  
  /// kinematic parameters of the fitted state
  TVectorT<double> kinparm_; //
  /// kinematic parameters of the unfitted state
  TVectorT<double> input_kinparm_; //
  
  /// coveriance matrix of the fitted state
  TMatrixDSym kinmatrix_; //
  /// coveriance matrix of the unfitted state
  TMatrixDSym input_kinmatrix_; //
  
  /// reference to the initial candidate
  reco::RecoChargedCandidateRef CandRef_; //
};

typedef std::vector<SelectedKinematicParticle> SelectedKinematicParticleCollection;

#endif
