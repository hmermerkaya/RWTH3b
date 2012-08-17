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
// $Id: SelectedKinematicDecay.h,v 1.21 2012/08/15 09:56:27 inugent Exp $
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
  enum FitSequence{ZeroAmbiguitySolution=0,PlusSolution,MinusSolution,NAmbiguity};

  SelectedKinematicDecay();
  SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRefOrig, std::vector<reco::TrackRef> &TrackTriplet, 
			 const reco::Vertex &primaryVertex,std::string primVtxReFitTag, unsigned int nTauPerVtx);

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
  void SetInitalKinematics(TVector3 tauFlghtDirNoCorr,std::vector<TLorentzVector> initalpions,TLorentzVector intial_a1_p4,
			   TVector3 tauFlghtDir,double initThetaGJ, double ThetaMax);
  void SetKinematicFitStatus(unsigned int ambiguity, const std::map<std::string, bool> discriminators);
  void SetKinematicFitProperties(unsigned int ambiguity,const SelectedKinematicParticleCollection particles, 
				 const int iterations, const int maxiterations, const float csum, 
				 const float mincsum, const int constraints, const int ndf, 
				 const float chi2);
  /// store quality discriminators that cannot directly be calculated from stored members only (e.g. conversion into reco::Vertex format would be needed). FIXME: replace this by  a dynamic calculation (depending on the decay mode)
  void SetQualityCriteria(unsigned int ambiguity, const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction); 
  void SetInitalGuess(unsigned int ambiguity,TLorentzVector &TauGuessLV,TLorentzVector &NuGuessLV);
  void SetKFSecondaryVertex(unsigned int ambiguity,reco::Vertex SecVtx);

  //////////////////////////////////////////////////////////////////////////
  //
  // Get Functions
  //
  unsigned int                        TauDecayMode(){return tauDecayMode_;}
  unsigned int                        NumberOfTauPerVtx()const{return nTauPerVtx_;}
  std::string                         PrimaryVertexReFitCollectionTag(){return primVtxReFitTag_;}

  // Inital variables
  const reco::PFTauRef              & PFTauRef() const { return PFTauRefOrig_;}
  const int                           sgnlConeTrkSize() const;                                  /// size of tracks in the signal cone of the initial PFTau candidate     
  const std::vector<reco::TrackRef> & InitalTrackTriplet()const{return initalTrackTriplet_;}
  const reco::Vertex                & InitalPrimaryVertex() const{return initalPrimVtx_;}
  std::vector<reco::Track>            InitalSecondaryVertexTracks(){return initalSecVtxTracks_;}
  reco::Vertex                        InitalSecondaryVertex(){return initalSecVtx_;};
  const reco::Vertex                & InitalPrimaryVertexReFit()const{return initalPrimaryVertexReFit_;}
  const reco::Vertex                & InitalPrimaryVertexReFitAndRotated()const{return initalPrimaryVertexReFitAndRotated_;}
  TLorentzVector                      Inital_a1_p4();
  double                              InitalThetaMax(){return initalThetaMax_;}
  double                              InitalThetaGJ(){return initalThetaGJ_;}
  TVector3                            InitalTauFlightDirectionReFitandRotatedPvtxToSvtx(){return initalTauFlghtDir_;}
  TVector3                            InitalTauFlightDirectionReFitPvtxToSvtx(){return initalTauFlghtDirNoCorr_;}
  TLorentzVector                      InitialTauGuess(unsigned int ambiguity);
  TLorentzVector                      InitalNeutrinoGuess(unsigned int ambiguity);
  std::vector<TLorentzVector>         InitalPions();
  const reco::Vertex                & PrimaryVertexReFitAndRotated()const{return initalPrimaryVertexReFitAndRotated_;}

  ///////////////////////////////////////////////////////////////////////
  // KF Info
  TLorentzVector                         Tau(unsigned int ambiguity=ZeroAmbiguitySolution);
  TLorentzVector                         Neutrino(unsigned int ambiguity=ZeroAmbiguitySolution);
  std::vector<TLorentzVector>            Pions(unsigned int ambiguity=ZeroAmbiguitySolution);
  TLorentzVector                         a1_p4(unsigned int ambiguity=ZeroAmbiguitySolution);
  reco::Vertex                           SecondaryVertex(unsigned int ambiguity=ZeroAmbiguitySolution);
      
  // KF variables
  const std::map<std::string,bool>            discriminators(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const SelectedKinematicParticleCollection   particles(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const SelectedKinematicParticle           * topParticle(unsigned int ambiguity=ZeroAmbiguitySolution) const;                       /// return the mother particle of this decay. ie the tau
  void                                        daughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity=ZeroAmbiguitySolution) const;        /// return all daughter particles 
  void                                        chargedDaughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity=ZeroAmbiguitySolution) const; /// return only all charged daughters particles 
  void                                        neutralDaughters(std::vector< SelectedKinematicParticle const * > & par,unsigned int ambiguity=ZeroAmbiguitySolution) const; /// return only all neutral daughter particles 
  void                                        modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par,unsigned int ambiguity=ZeroAmbiguitySolution); ///DO NOT USE after reading from event stream!
  
  const int    iterations(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const int    maxiterations(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const float  chi2(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const float  constraints(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const float  ndf(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const float  csum(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const float  mincsum(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const double vtxSignPVRotSV(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const double vtxSignPVRotPVRed(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const double a1Mass(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  const double chi2prob(unsigned int ambiguity=ZeroAmbiguitySolution) const;       
  const double energyTFraction(unsigned int ambiguity=ZeroAmbiguitySolution) const;
  


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
  TVector3                     initalTauFlghtDirNoCorr_;
  std::vector<TLorentzVector>  initalTauGuess_;
  std::vector<TLorentzVector>  initalNuGuess_;
  std::vector<TLorentzVector>  initalpions_;
  TLorentzVector               intial_a1_p4_;

  // KF variables
  std::vector<reco::Vertex>                 SecVtx_;             /// secondary vertex from kinematic fit
  SelectedKinematicParticleCollection       particles_;          /// collection of kinematic particles assigned to this decay
  std::vector<int>                          iterations_;         /// fit parameter: iterations until convergence
  std::vector<int>                          maxiterations_;      /// fit parameter: maximal allowed iterations
  std::vector<float>                        csum_;               /// fit parameter: sum of constraints after last iteration
  std::vector<float>                        mincsum_;            /// fit parameter: minimal sum of constraints. fall below for convergence.
  std::vector<float>                        chi2_;               /// fit parameter: chi2 of the fit 
  std::vector<float>                        constraints_;        /// fit parameter: number of constraints applied to the fit. (WARNING: This values is called ndf in RecoVertex/KinematicFit)
  std::vector<float>                        ndf_;                /// fit parameter: real ndf depending on the decay mode. (WARNING: This is NOT the ndf() in RecoVertex/KinematicFit)
  std::vector<std::map<std::string, bool> > discriminators_;     /// official pftau discriminators (ensure same size and order)
  
  /// quality criterion (may depend on decay mode)
  std::vector<double>                       vtxSignPVRotSV_;     /// vertex significance between the rotated primary vertex and the secondary vertex of the tau decay 
  std::vector<double>                       vtxSignPVRotPVRed_;  /// vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
  std::vector<double>                       a1Mass_;             /// mass of the a1 system
  std::vector<double>                       energyTFraction_;    /// transversal energy fraction between the intial PFTau and the final kinematic tau

};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
