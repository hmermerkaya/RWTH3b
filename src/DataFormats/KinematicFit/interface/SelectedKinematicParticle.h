#ifndef DataFormats_KinematicFit_SelectedKinematicParticle_h
#define DataFormats_KinematicFit_SelectedKinematicParticle_h
// -*- C++ -*-
//
// Package:    SelectedKinematicParticle
// Class:      SelectedKinematicParticle
// 
/**
 
 Description: saves results from kinematically refit particles
 
 Implementation:
 <Notes on implementation>
 */
//
// $Id: SelectedKinematicParticle.h,v 1.4 2010/01/20 15:32:55 perchall Exp $
//
//


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
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


class SelectedKinematicParticle {
public:
	SelectedKinematicParticle();
	SelectedKinematicParticle(const TVectorT<double> kinparm_, const TMatrixDSym kinmatrix_, const int charge_, const std::string name_, const float chi2_, const float ndf_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_, const int status_);
	SelectedKinematicParticle(const TVectorT<double> kinparm_, const TMatrixDSym kinmatrix_, const TVectorT<double> input_kinparm_, const TMatrixDSym input_kinmatrix_, const int charge_, const std::string name_, const float chi2_, const float ndf_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_, const int status_);
    SelectedKinematicParticle(const RefCountedKinematicParticle kinparticle_, const std::string name_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_, const int status_);

    int status() const;
    int matched() const;
    int iterations() const;
    int maxiterations() const;
    int ambiguity() const;

    int charge() const;
	float chi2() const;
	float ndf() const;
    float csum() const;
    float mincsum() const;

    std::string name() const;
    
	TVectorT<double> parameters() const;
	TVectorT<double> input_parameters() const;

	TMatrixDSym matrix() const;
	TMatrixDSym input_matrix() const;
    
	reco::RecoChargedCandidateRef candRef() const;
	void setCandRef(const reco::RecoChargedCandidateRef parm);
    
    TLorentzVector p4() const;
	TVector3 vertex() const;
    
	void setMatched(const int parm);
	void setInitialTauState(const TLorentzVector & tau, const reco::Vertex & primVtx);//initial tau state consists of used primVtx (including errors) and the visible tau parameters
	
private:
    int status__; //
    int matched__; //
    int iterations__; //
    int maxiterations__; //
    int ambiguity__; //
    
    int charge__; //
    float chi2__; // 
	float ndf__; //
    float csum__; //
    float mincsum__; //
    
    std::string name__; //
    
	TVectorT<double> kinparm__; //
	TVectorT<double> input_kinparm__; //

	TMatrixDSym kinmatrix__; //
	TMatrixDSym input_kinmatrix__; //
    
    reco::RecoChargedCandidateRef CandRef__; //
    
    TVectorT<double> convertVector( const AlgebraicVector7 vector );
    TMatrixDSym convertMatrix( const AlgebraicSymMatrix77 matrix );
};

typedef std::vector<SelectedKinematicParticle> SelectedKinematicParticleCollection;

#endif
