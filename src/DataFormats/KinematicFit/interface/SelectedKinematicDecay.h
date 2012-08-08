#ifndef DataFormats_KinematicFit_SelectedKinematicDecay_h
#define DataFormats_KinematicFit_SelectedKinematicDecay_h

// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      SelectedKinematicDecay
// 
/**
 * This data format combines all objects of type SelectedKinematicParticle assigned to one tau decay.
 * WARNING: the current implementation of quality cuts may assume a certain decay mode!!!
 *
 * @author Lars Perchalla, Philip Sauerland
 * @date 2010
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Jan  21 17:29:43 CEST 2010
// $Id: SelectedKinematicDecay.h,v 1.19 2012/08/01 09:22:13 inugent Exp $
//
//

#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"
#include <string>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class SelectedKinematicDecay {
 public:
  enum TauType{Undefined=0,ThreePion,ThreePionAndPi0};
  enum FitSequence{AmbiguitySolution=0,PlusSolution,MinusSolution};

  SelectedKinematicDecay();
  SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRefOrig, std::vector<reco::TrackRef> &TrackTriplet, 
			 const reco::Vertex &primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx);
  SelectedKinematicDecay(const SelectedKinematicParticleCollection & particles, 
			 const int iterations, const int maxiterations, const float csum, 
			 const float mincsum, const int constraints, const int ndf, const float chi2, 
			 const reco::PFTauRef & tauRef, const std::map<std::string, bool> & discriminators);

  virtual ~SelectedKinematicDecay();

  //////////////////////////////////////////////////////////////////////////
  //
  // Set Functions
  // 
  void SetTauDecayMode(unsigned int tauDecayMode){tauDecayMode_=tauDecayMode;}
  void SetInitialProperties(unsigned int tauDecayMode, const reco::PFTauRef tauRefOrig, std::vector<reco::TrackRef> TrackTriplet, 
			    const reco::Vertex primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx);
  void SetInitalVertexProperties(reco::Vertex primaryVertexReFit,reco::Vertex primaryVertexReFitAndRotated,
			      std::vector<reco::TransientTrack> secVtxTracks, TransientVertex secVtx);
  void SetInitalKinematics(TVector3 tauFlghtDir,std::vector<TLorentzVector> initalpions,TLorentzVector intial_a1_p4_,double initThetaGJ, double ThetaMax);
  void SetKinematicFitProperties(const SelectedKinematicParticleCollection particles, 
				 const int iterations, const int maxiterations, const float csum, 
				 const float mincsum, const int constraints, const int ndf, 
				 const float chi2, const std::map<std::string, bool> discriminators);
  /// store quality discriminators that cannot directly be calculated from stored members only (e.g. conversion into reco::Vertex format would be needed). FIXME: replace this by  a dynamic calculation (depending on the decay mode)
  void SetMissingQualityCriteria(const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction); 
  void SetInitalGuess(std::vector<TLorentzVector> &TauGuessLV,std::vector<TLorentzVector> &NuGuessLV);
  void SetKFSecondaryVertex(reco::Vertex SecVtx_);

  //////////////////////////////////////////////////////////////////////////
  //
  // Get Functions
  //
  unsigned int                        TauDecayMode(){return tauDecayMode_;}
  unsigned int                        NumberOfTauPerVtx()const{return nTauPerVtx_;}
  std::string                         PrimaryVertexReFitCollectionTag(){return primVtxReFitTag_;}

  // Inital variables
  const reco::PFTauRef              & PFTauRef() const { return PFTauRefOrig_;}
  const std::vector<reco::TrackRef> & InitalTrackTriplet()const{return initalTrackTriplet_;}
  const reco::Vertex                & InitalPrimaryVertex() const{return initalPrimVtx_;}
  std::vector<reco::Track>            InitalSecondaryVertexTracks(){return initalSecVtxTracks_;}
  reco::Vertex                        InitalSecondaryVertex(){return initalSecVtx_;};
  const reco::Vertex                & InitalPrimaryVertexReFit()const{return initalPrimaryVertexReFit_;}
  const reco::Vertex                & InitalPrimaryVertexReFitAndRotated()const{return initalPrimaryVertexReFitAndRotated_;}
  TLorentzVector                      Inital_a1_p4();
  double                              InitalThetaMax(){return initalThetaMax_;}
  double                              InitalThetaGJ(){return initalThetaGJ_;}
  TVector3                            InitalFlightDirection(){return initalTauFlghtDir_;}
  TLorentzVector                      InitialTauGuess(unsigned int ambiguity);
  TLorentzVector                      InitalNeutrinoGuess(unsigned int ambiguity);
  std::vector<TLorentzVector>         InitalPions();

  // KF Info
  TLorentzVector                         Tau(unsigned int ambiguity);
  TLorentzVector                         Neutrino(unsigned int ambiguity);
  std::vector<TLorentzVector>            Pions(unsigned int ambiguity);
  const reco::Vertex                   & PrimaryVertexReFitAndRotated()const{return initalPrimaryVertexReFitAndRotated_;}
  reco::Vertex                           SecondaryVertex(unsigned int ambiguity){return SecVtx_;}
      
  // KF variables
  const std::map<std::string, bool>         & discriminators()const{ return discriminators_; }
  const SelectedKinematicParticleCollection & particles() const { return particles_; } /// return all particles assigned to this decay including the mother.
  const SelectedKinematicParticle*            topParticle() const;                                /// return the mother particle of this decay. this is the decaying particle itself.
  void                                        daughters(std::vector< SelectedKinematicParticle const * > & par) const;        /// return all particles assigned to this decay w/o the mother
  void                                        chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const; /// return only all charged particles assigned to this decay w/o the mother
  void                                        neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const; /// return only all neutral particles assigned to this decay w/o the mother
  void                                        modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par); ///DO NOT USE after reading from event stream!
  
  const int    iterations() const{return iterations_;}
  const int    maxiterations() const{return maxiterations_;}
  const float  chi2() const{return chi2_;}
  const float  constraints() const{return constraints_;}
  const float  ndf() const{return ndf_;}
  const float  csum() const{return csum_;}
  const float  mincsum() const{return mincsum_;}
  const double vtxSignPVRotSV() const { return vtxSignPVRotSV_; }        /// quality criterion: vertex significance between the rotated primary vertex and the secondary vertex of the tau decay
  const double vtxSignPVRotPVRed() const { return vtxSignPVRotPVRed_; }  /// quality criterion: vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
  const double a1Mass() const { return a1Mass_;}                         /// quality criterion: mass of the a1 system
  const double chi2prob() const;                                         /// quality criterion: chi2 probability of the fit
  const int    sgnlConeTrkSize() const;                                  /// size of tracks in the signal cone of the initial PFTau candidate
  const double energyTFraction() const { return energyTFraction_; }      /// quality criterion: transversal energy fraction between the intial PFTau and the final kinematic tau
  


 private:
  // internal variables
  unsigned int                 tauDecayMode_;
  reco::PFTauRef               PFTauRefOrig_;
  std::vector<reco::TrackRef>  initalTrackTriplet_;
  reco::Vertex                 initalPrimVtx_;
  std::string                  primVtxReFitTag_;
  unsigned int                 nTauPerVtx_;
  std::vector<reco::Track>     initalSecVtxTracks_; 
  reco::Vertex                 initalSecVtx_;
  reco::Vertex                 initalPrimaryVertexReFit_;
  reco::Vertex                 initalPrimaryVertexReFitAndRotated_;
  double                       initalThetaMax_;
  double                       initalThetaGJ_;
  TVector3                     initalTauFlghtDir_;
  std::vector<TLorentzVector>  initalTauGuess_;
  std::vector<TLorentzVector>  initalNuGuess_;
  std::vector<TLorentzVector>  initalpions_;
  TLorentzVector               intial_a1_p4_;

  // KF variables
  reco::Vertex                        SecVtx_;    /// secondary vertex from kinematic fit
  SelectedKinematicParticleCollection particles_; /// collection of kinematic particles assigned to this decay
  int iterations_;                                /// fit parameter: iterations until convergence
  int maxiterations_;                             /// fit parameter: maximal allowed iterations
  float csum_;                                    /// fit parameter: sum of constraints after last iteration
  float mincsum_;                                 /// fit parameter: minimal sum of constraints. fall below for convergence.
  float chi2_;                                    /// fit parameter: chi2 of the fit 
  float constraints_;                             /// fit parameter: number of constraints applied to the fit. (WARNING: This values is called ndf in RecoVertex/KinematicFit)
  float ndf_;                                     /// fit parameter: real ndf depending on the decay mode. (WARNING: This is NOT the ndf() in RecoVertex/KinematicFit)
  std::map<std::string, bool> discriminators_;    /// official pftau discriminators (ensure same size and order)
  
  /// quality criterion (may depend on decay mode)
  double vtxSignPVRotSV_;                         /// vertex significance between the rotated primary vertex and the secondary vertex of the tau decay 
                                                  /// (the tau carries rotated and reduced primVtx as initial vtx, but conversion into reco::Vertex needed)
  double vtxSignPVRotPVRed_;                      /// vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
  double a1Mass_;                                 /// mass of the a1 system
  double energyTFraction_;                        /// transversal energy fraction between the intial PFTau and the final kinematic tau

};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
