// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      Tau_JAKID_Filter
// 
/**\class Tau_JAKID_Filter Tau_JAKID_Filter.cc RecoTauTag/KinematicTau/src/Tau_JAKID_Filter.cc
 
 Description: Filters on JAK ID
Added by: Ian M. Nugent
Aug 19 2012  
 */

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


class Tau_JAKID_Filter : public edm::EDFilter {
public:
  explicit Tau_JAKID_Filter(const edm::ParameterSet&);
  ~Tau_JAKID_Filter();
	
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual bool isTruthTauInAcceptance(const reco::GenParticle &cand);

  std::vector<int> JAKID_;
  std::vector<int> nprongs_;
  edm::InputTag gensrc_;
  double TauPtMin_,TauEtaMax_; 
  int cnt_,cntFound_; 


};
