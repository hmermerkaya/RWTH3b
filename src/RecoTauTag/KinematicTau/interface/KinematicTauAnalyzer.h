/**
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010
Modified by Ian M. Nugent 
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
  std::vector<int> JAKID_;
  int cnt_;
  std::vector<int> cntFound_;

  /////////////////////////////////
  //
  // Valdiation Histograms
  //
  // Check constraints
  DQMStore *dbe;
  std::vector<MonitorElement*> 
    nEvt, 
    TauMass, 
    PionMass, 
    NuMass, 
    dTauMass, 
    dPionMass, 
    dNuMass,
    VtxXChange, 
    VtxYChange, 
    VtxZChange, 
    SecVtxXChange, 
    SecVtxYChange, 
    SecVtxZChange,
    TauPhiChange, 
    TauThetaChange, 
    TauEChange, 
    PionPhiChange, 
    PionThetaChange, 
    PionEChange,
    NuPhiChange, 
    NuThetaChange, 
    NuEChange, 
    JAKID, 
    JAKIDall, 
    JAKIDeff, 
    Truth_TauMatched,
    vtxSignPVRotSV, 
    vtxSignPVRotPVRed, 
    a1Mass, 
    energyTFraction, 
    iterations, 
    maxiterations,
    chi2, 
    constraints, 
    ndf, 
    csum, 
    mincsum, 
    chi2prob, 
    TauFlightDir, 
    TauFlightDirInitial,
    GFAngleInitial,
    GFAngle;

  std::vector<std::vector<MonitorElement*> >  
    Truth_TauMatch_dPhi, 
    Truth_TauMatch_dTheta, 
    Truth_TauMatch_dE, 
    Truth_PionMatch_dPhi, 
    Truth_PionMatch_dTheta, 
    Truth_PionMatch_dE,
    Truth_NuMatch_dPhi,
    Truth_NuMatch_dTheta,
    Truth_NuMatch_dE,
    Truth_A1Match_dPhi, 
    Truth_A1Match_dTheta, 
    Truth_A1Match_dE, 
    Truth_A1Match_dPt, 
    Truth_A1Match_M,
    TruthVtxX, 
    TruthVtxY, 
    TruthVtxZ, 
    TruthSecVtxX, 
    TruthSecVtxY, 
    TruthSecVtxZ, 
    PullVtxX, 
    PullVtxY, 
    PullVtxZ, 
    PullSecVtxX, 
    PullSecVtxY, 
    PullSecVtxZ,
    TruthPVtxSig, 
    TruthSecVtxSig, 
    TruthTauFlightDir, 
    TruthTauFlightDirCheck,
    Truth_TauMatch_dPt,
    Truth_TauMatch_dPz, 
    Truth_TauMatch_dPtvsL, 
    Truth_TauMatch_dEvsL,
    Truth_TauMatch_dphivsL,
    Truth_TauMatch_dthetavsL, 
    Truth_TauMatch_dGFAnglevsL, 
    Truth_TauMatch_dGFAngle, 
    Truth_A1Match_dPtvslength, 
    Truth_A1Match_dphivslength, 
    Truth_A1Match_dthetavslength,
    Truth_PionMatch_dPtvslength, 
    Truth_PionMatch_dphivslength, 
    Truth_PionMatch_dthetavslength,
    Truth_TauMatch_TrueMaxGFAngle,
    Truth_TauMatch_TrueGFAngle,
    Truth_TauMatch_TrueGFAngleoverMaxGFAngle,
    Truth_TauMatch_dGFAnglevsTrackchi2,
    Truth_TauMatch_dGFAnglevsTrackQuality,
    Truth_TauMatch_dGFAnglevsMaxTrackPt,
    Truth_TauMatch_dGFAnglevsMinTrackPt,
    Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt,
    Truth_TauMatch_dGFAnglevsHitFrac,
    PullTauPx,
    PullTauPy,
    PullTauPz;



  std::map<unsigned int,unsigned int> JAKIDtoIndex;
 

};

#endif
