/**
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */

#ifndef KinematicTauAnalyzer_h
#define KinematicTauAnalyzer_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <map>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


class KinematicTauAnalyzer : public edm::EDAnalyzer {
 public:
  explicit KinematicTauAnalyzer(const edm::ParameterSet&);
  ~KinematicTauAnalyzer();
  
private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual bool doJAKID(unsigned int i);
  virtual bool isTruthTauInAcceptance(const reco::GenParticle &cand);
  
  std::vector<std::string> discriminators_;
  edm::InputTag KinematicFitTauTag_,gensrc_,GenEventInfo_;
  float TauMatchingDR_,TauPtMin_,TauEtaMax_,tau_pdgid;
  int tau_pdgid_;
  unsigned int NJAKID_;

  /////////////////////////////////
  //
  // Valdiation Histograms
  //
  // Check constraints
  DQMStore *dbe;
  std::vector<MonitorElement*> nEvt, TauMass, PionMass, NuMass, dTauMass, dPionMass, dNuMass,
    VtxXChange, VtxYChange, VtxZChange, SecVtxXChange, SecVtxYChange, SecVtxZChange,
    TauPhiChange, TauThetaChange, TauEChange, PionPhiChange, PionThetaChange, PionEChange,
    NuPhiChange, NuThetaChange, NuEChange, JAKID, JAKIDall, JAKIDeff, TauMatched;
  std::vector<std::vector<MonitorElement*> >  TauMatch_dphi, TauMatch_dtheta, TauMatch_e, 
    PionMatch_dphi, PionMatch_dtheta, PionMatch_e, NuMatch_dphi, NuMatch_dtheta, NuMatch_e;
  std::map<unsigned int,unsigned int> JAKIDtoIndex;
 

};

#endif
