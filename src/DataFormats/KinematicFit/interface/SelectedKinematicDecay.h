#ifndef DataFormats_KinematicFit_SelectedKinematicDecay_h
#define DataFormats_KinematicFit_SelectedKinematicDecay_h
// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      SelectedKinematicDecay
// 
/**
 * This data format combines all objects of type SelectedKinematicParticle assigned to one tau decay.
 * @author Lars Perchalla, Philip Sauerland
 * @date 2010
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Jan  21 17:29:43 CEST 2010
// $Id: SelectedKinematicDecay.h,v 1.14 2011/12/08 16:32:28 sauerlan Exp $
//
//

#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"


class SelectedKinematicDecay {
public:
	SelectedKinematicDecay();
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles);
	SelectedKinematicDecay(SelectedKinematicParticleCollection particles, const int & signalPFChargedHadrCands, const int & signalPFNeutrHadrCands, 
			       const double & fraction, const double & RefitMass, const double & Chi2, const double & ModPV_PV_significance,const double & SV_PV_significance, const std::map<std::string, bool> & discriminators);

    /**
	 return the mother particle of this decay. this is the decaying particle itself.
	 */
	const SelectedKinematicParticle* topParticle() const;
    /**
	 return all particles assigned to this decay including the mother.
	 */
    void particles(std::vector< SelectedKinematicParticle const * > & par) const;
	/**
	 return all particles assigned to this decay w/o the mother
	 */	
    void daughters(std::vector< SelectedKinematicParticle const * > & par) const;
	/**
	 return only all charged particles assigned to this decay w/o the mother
	 */	
    void chargedDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
	/**
	 return only all neutral particles assigned to this decay w/o the mother
	 */	
    void neutralDaughters(std::vector< SelectedKinematicParticle const * > & par) const;
	/**
	 return the number of charged candidates found within the signal cone of the initial PFlow candidate
	 */	
    int signalPFChargedHadrCands() const ;
	/**
	 return the number of neutral candidates found within the signal cone of the initial PFlow candidate
	 */	
    int signalPFNeutrHadrCands() const ;
	/**
	 return a map of discriminators of the initial PFlow candidate
	 */	
    double TauEnergyFraction() const ;
	/**
	 return the energy ratio
	 */	
    double RefitVisibleMass() const ;
	/**
	 return reffited visible mass
	 */	
    double Chi2() const ;
	/**
	 return chi2 of the fit
	 */	
    double ModPV_PV_significance() const ;
	/**
	 return a map of discriminators of the initial PFlow candidate
	 */	
    double PV_SV_significance() const ;
	/**
	 return a map of discriminators of the initial PFlow candidate
	 */	

    std::map<std::string, bool> discriminators() const;
	/**
	 DO NOT USE after reading from event stream!
	 */
    void modifiableChargedDaughters(std::vector< SelectedKinematicParticle * > & par);
 	
	
	void setPFTauRef(const std::vector<reco::PFTauRef> & value);
	void addPFTauRef(const reco::PFTauRef & value);
	/**
	 index of the initial PFTau from which this decay was created. multiple indeces possible
	 */
	std::vector<reco::PFTauRef> PFTauRef() const;

private:
    SelectedKinematicParticleCollection particles_;
    ///size of tracks in the signal cone of the PFTau candidate
	int signalPFChargedHadrCands_;
    ///size of neutral candidates in the signal cone of the PFTau candidate
	int signalPFNeutrHadrCands_;
    ///official pftau discriminators (ensure same size and order)
	std::map<std::string, bool> discriminators_;
    ///store the index of the initial PFTau from which this decay was created. multiple indeces possible
	std::vector<reco::PFTauRef> PFTauRef_;

   ///Quality criteria
	double fraction_;
        double RefitMass_;
        double Chi2_;
        double ModPV_PV_significance_;
	double SV_PV_significance_;
};

typedef std::vector<SelectedKinematicDecay> SelectedKinematicDecayCollection;

#endif
