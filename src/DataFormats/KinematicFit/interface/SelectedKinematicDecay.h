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
// $Id: SelectedKinematicDecay.h,v 1.15 2012/01/20 14:17:44 cherepan Exp $
//
//

#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

class SelectedKinematicDecay {
public:
	SelectedKinematicDecay();
	SelectedKinematicDecay(const SelectedKinematicParticleCollection & particles, const int iterations, const int maxiterations, const float csum, const float mincsum, const int constraints, const int ndf, const float chi2, const reco::PFTauRef & tauRef, const std::map<std::string, bool> & discriminators/*, const reco::VertexRef & primaryVertexRef*/);

	const reco::PFTauRef & PFTauRef() const { return PFTauRef_; }
    const std::map<std::string, bool> & discriminators() const { return discriminators_; }
//	const reco::VertexRef & primaryVertexRef() const { return primaryVertexRef_; }

	/// return all particles assigned to this decay including the mother.
    const SelectedKinematicParticleCollection & particles() const { return particles_; }
    /// return the mother particle of this decay. this is the decaying particle itself.
	const SelectedKinematicParticle* topParticle() const;
	/// return all particles assigned to this decay w/o the mother
    void daughters(std::vector< SelectedKinematicParticle const * > & par) const;
	/// return only all charged particles assigned to this decay w/o the mother
    void chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
	/// return only all neutral particles assigned to this decay w/o the mother
    void neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
    
	/**
	 DO NOT USE after reading from event stream!
	 */
    void modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par);
    
    const int iterations() const;
    const int maxiterations() const;
	const float chi2() const;
    const float constraints() const;
	const float ndf() const;
    const float csum() const;
    const float mincsum() const;
    
    /// store quality discriminators that cannot directly be calculated from stored members only (e.g. conversion into reco::Vertex format would be needed). FIXME: replace this by  a dynamic calculation (depending on the decay mode)
    void setMissingQualityCriteria(const double & vtxSignPVRotSV, const double & vtxSignPVRotPVRed, const double & a1Mass, const double & energyTFraction);
    /// quality criterion: vertex significance between the rotated primary vertex and the secondary vertex of the tau decay
    const double vtxSignPVRotSV() const { return vtxSignPVRotSV_; }
    /// quality criterion: vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
    const double vtxSignPVRotPVRed() const { return vtxSignPVRotPVRed_; }
    /// quality criterion: mass of the a1 system
    const double a1Mass() const { return a1Mass_; }
    /// quality criterion: chi2 probability of the fit
    const double chi2prob() const;
    /// size of tracks in the signal cone of the initial PFTau candidate
    const int sgnlConeTrkSize() const;
    /// quality criterion: transversal energy fraction between the intial PFTau and the final kinematic tau
    const double energyTFraction() const { return energyTFraction_; }
    
private:
    /// collection of kinematic particles assigned to this decay
    SelectedKinematicParticleCollection particles_;
    /// fit parameter: iterations until convergence
    int iterations_; //
    /// fit parameter: maximal allowed iterations
    int maxiterations_; //
    /// fit parameter: sum of constraints after last iteration
    float csum_; //
    /// fit parameter: minimal sum of constraints. fall below for convergence.
    float mincsum_; //
    /// fit parameter: chi2 of the fit
    float chi2_; // 
    /// fit parameter: number of constraints applied to the fit. (WARNING: This values is called ndf in RecoVertex/KinematicFit)
	float constraints_; //
    /// fit parameter: real ndf depending on the decay mode. (WARNING: This is NOT the ndf() in RecoVertex/KinematicFit)
	float ndf_; //
    
    /// reference of the initial PFTau from which this decay was created
	reco::PFTauRef PFTauRef_;
    /// official pftau discriminators (ensure same size and order)
	std::map<std::string, bool> discriminators_;
//    /// reference of the unrotated primary vertex
//    reco::VertexRef primaryVertexRef_;
    
    /// quality criterion (may depend on decay mode): vertex significance between the rotated primary vertex and the secondary vertex of the tau decay (the tau carries rotated and reduced primVtx as initial vtx, but conversion into reco::Vertex needed)
    double vtxSignPVRotSV_;
    /// quality criterion (may depend on decay mode): vertex significance of the primary vertex rotation w.r.t. the initial primary vertex (already w/o tracks assigned to the tau decay)
    double vtxSignPVRotPVRed_;
    /// quality criterion (may depend on decay mode): mass of the a1 system
    double a1Mass_;
    /// quality criterion (may depend on decay mode): transversal energy fraction between the intial PFTau and the final kinematic tau
    double energyTFraction_;

};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
