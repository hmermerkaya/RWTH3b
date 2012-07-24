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
// $Id: SelectedKinematicDecay.h,v 1.16 2012/01/25 12:58:11 perchall Exp $
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

class SelectedKinematicDecay {
 public:
  enum TauType{Undefined=0,ThreePion,ThreePionAndPi0};
  enum FitSequence{AmbiguitySolution=0,PlusSolution,MinusSolution};

  SelectedKinematicDecay();
  SelectedKinematicDecay(unsigned int tauDecayMode, const reco::PFTauRef &tauRef, reco::TrackRefVector &TrackTriplet, 
			 const reco::VertexRef &primaryVertexRef,std::string primVtxReFitTag, unsigned int nTauPerVtx,
			 std::vector<reco::TransientTrack> &secVtxTracks, TransientVertex &secVtx);
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
  void SetInitialProperties(unsigned int tauDecayMode, const reco::PFTauRef tauRef,reco::TrackRefVector TrackTriplet, 
			    const reco::VertexRef primaryVertexRef,std::string primVtxReFitTag, unsigned int nTauPerVtx,
			    std::vector<reco::TransientTrack> secVtxTracks, TransientVertex secVtx);
  void SetPrimaryVertexReFit(reco::VertexRef primaryVertexReFit);
  void SetPrimaryVertexReFitAndRotated(reco::Vertex primaryVertexReFitAndRotated);
  void SetKinematicFitProperties(const SelectedKinematicParticleCollection particles, 
				 const int iterations, const int maxiterations, const float csum, 
				 const float mincsum, const int constraints, const int ndf, 
				 const float chi2, const std::map<std::string, bool> discriminators);
  /// store quality discriminators that cannot directly be calculated from stored members only (e.g. conversion into reco::Vertex format would be needed). FIXME: replace this by  a dynamic calculation (depending on the decay mode)
  void setMissingQualityCriteria(const double vtxSignPVRotSV, const double vtxSignPVRotPVRed, const double a1Mass, const double energyTFraction); 

  //////////////////////////////////////////////////////////////////////////
  //
  // Get Functions
  //
  unsigned int TauDecayMode(){return tauDecayMode_;}
  const reco::PFTauRef & PFTauRef() const { return PFTauRef_;}
  const reco::TrackRefVector & TrackTriplet()const{return TrackTriplet_;}
  const reco::VertexRef & PrimaryVertex() const{return primVtx_;}
  std::string PrimaryVertexReFitCollectionTag(){return primVtxReFitTag_;}
  unsigned int NumberOfTauPerVtx()const{return nTauPerVtx_;}
  std::vector<reco::Track> SecondaryVertexTracks(){return secVtxTracks_;}
  reco::Vertex  SecondaryVertex(){return secVtx_;};

  const reco::VertexRef & PrimaryVertexReFit()const{return primaryVertexReFit_;}
  const reco::Vertex & PrimaryVertexReFitAndRotated()const{return primaryVertexReFitAndRotated_;}
  const std::map<std::string, bool> & discriminators() const { return discriminators_; }
    
  const SelectedKinematicParticleCollection & particles() const { return particles_; } /// return all particles assigned to this decay including the mother.
  const SelectedKinematicParticle* topParticle() const;                                /// return the mother particle of this decay. this is the decaying particle itself.
  void daughters(std::vector< SelectedKinematicParticle const * > & par) const;        /// return all particles assigned to this decay w/o the mother
  void chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const; /// return only all charged particles assigned to this decay w/o the mother
  void neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const; /// return only all neutral particles assigned to this decay w/o the mother
  
  void modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par); ///DO NOT USE after reading from event stream!
  
  const int iterations() const;
  const int maxiterations() const;
  const float chi2() const;
  const float constraints() const;
  const float ndf() const;
  const float csum() const;
  const float mincsum() const;
  const double vtxSignPVRotSV() const { return vtxSignPVRotSV_; }        /// quality criterion: vertex significance between the rotated primary vertex and the secondary vertex of the tau decay
  const double vtxSignPVRotPVRed() const { return vtxSignPVRotPVRed_; }  /// quality criterion: vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
  const double a1Mass() const { return a1Mass_; }                        /// quality criterion: mass of the a1 system
  const double chi2prob() const;                                         /// quality criterion: chi2 probability of the fit
  const int sgnlConeTrkSize() const;                                     /// size of tracks in the signal cone of the initial PFTau candidate
  const double energyTFraction() const { return energyTFraction_; }      /// quality criterion: transversal energy fraction between the intial PFTau and the final kinematic tau
  
 private:
  // internal variables
  unsigned int             tauDecayMode_;
  reco::PFTauRef           PFTauRef_;
  reco::TrackRefVector     TrackTriplet_;
  reco::VertexRef          primVtx_;
  std::string              primVtxReFitTag_;
  unsigned int             nTauPerVtx_;
  std::vector<reco::Track> secVtxTracks_; 
  reco::Vertex             secVtx_;
  reco::VertexRef          primaryVertexReFit_;
  reco::Vertex             primaryVertexReFitAndRotated_;

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
                                                  ///(the tau carries rotated and reduced primVtx as initial vtx, but conversion into reco::Vertex needed)
  double vtxSignPVRotPVRed_;                      /// vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
  double a1Mass_;                                 /// mass of the a1 system
  double energyTFraction_;                        /// transversal energy fraction between the intial PFTau and the final kinematic tau

};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
