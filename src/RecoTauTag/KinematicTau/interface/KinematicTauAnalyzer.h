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
#include "TMatrixT.h"

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
  std::string tauType_;
  float TauMatchingDR_,TauPtMin_,TauEtaMax_,tau_pdgid;
  int tau_pdgid_;
  unsigned int NJAKID_;
  std::vector<int> JAKID_;
  int cnt_;
  std::vector<int> cntFound_;
  bool doFakeRate_, doDQM_;

  /////////////////////////////////
  //
  // Valdiation Histograms
  //
  // Check constraints

  DQMStore *dbe;
  std::vector<MonitorElement*> FakeRate_eta, FakeRate_pt,FakeRate_eta_All, FakeRate_pt_All, FakeRate_eta_Eff, FakeRate_pt_Eff, 
    FakeRate_pt_isFit, FakeRate_eta_isFit, FakeRate_pt_isFit_Eff, FakeRate_eta_isFit_Eff;

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
    TauPhi, 
    TauTheta, 
    TauE,
    TauPt,
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
    GFAngle_TauRF,
    GFAngle,
    chi2Vtx,
    ndfVtx,
    chi2probVtx,
    chi2SVtx,
    ndfSVtx,
    chi2probSVtx,
    FlightLength,
    FlightLengthSig,
    PionDr,
    PionDrHPS,
    PionDrHPSwithFit;

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
    TruthSecVtxXp,
    TruthSecVtxYp,
    TruthSecVtxZp,
    PullVtxX, 
    PullVtxY, 
    PullVtxZ, 
    PullSecVtxX, 
    PullSecVtxY, 
    PullSecVtxZ,
    PullSecVtxXp,
    PullSecVtxYp,
    PullSecVtxZp,
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
    Truth_TauMatch_reldPtvsL, 
    Truth_TauMatch_dGFAnglevsL, 
    Truth_TauMatch_dGFAngle,
    Truth_A1Match_dPtvslength, 
    Truth_A1Match_dphivslength, 
    Truth_A1Match_dthetavslength,
    Truth_A1Match_reldPtvslength,
    Truth_PionMatch_dPtvslength, 
    Truth_PionMatch_dphivslength, 
    Truth_PionMatch_dthetavslength,
    Truth_PionMatch_reldPtvslength,
    Truth_TauMatch_TrueMaxGFAngle,
    Truth_TauMatch_TrueGFAngle,
    Truth_TauMatch_TrueGFAngle_TauRF,
    Truth_TauMatch_GFAngleRecovsTruth_TauRF,
    Truth_TauMatch_TrueGFAngleoverMaxGFAngle,
    Truth_TauMatch_dGFAnglevsTrackchi2,
    Truth_TauMatch_dGFAnglevsTrackQuality,
    Truth_TauMatch_dGFAnglevsMaxTrackPt,
    Truth_TauMatch_dGFAnglevsMinTrackPt,
    Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt,
    Truth_TauMatch_dGFAnglevsHitFrac,
    Truth_TauMatch_dGFAnglevsVertexchi2,
    PullTauPx,
    PullTauPy,
    PullTauPz,
    Truth_TauMatch_dPtvsTheta,
    Truth_TauMatch_dEvsTheta,
    Truth_TauMatch_dthetavsTheta,
    Truth_TauMatch_reldPtvsTheta,
    Truth_TauMatch_dphivsTheta,
    Truth_TauMatch_dPtvsEta,
    Truth_TauMatch_dEvsEta,
    Truth_TauMatch_dthetavsEta,
    Truth_TauMatch_reldPtvsEta,
    Truth_TauMatch_dphivsEta,
    Truth_TauMatch_dPtvsPt,
    Truth_TauMatch_dEvsPt,
    Truth_TauMatch_dphivsPt,
    Truth_TauMatch_dthetavsPt,
    Truth_TauMatch_reldPtvsPt,
    Truth_TauMatch_dPtvsPtandEtaCut,
    Truth_TauMatch_dEvsPtandEtaCut,
    Truth_TauMatch_dphivsPtandEtaCut,
    Truth_TauMatch_dthetavsPtandEtaCut,
    Truth_TauMatch_reldPtvsPtandEtaCut,
    Truth_TauMatch_dPtvsLandEtaCut,
    Truth_TauMatch_dEvsLandEtaCut,
    Truth_TauMatch_dphivsLandEtaCut,
    Truth_TauMatch_dthetavsLandEtaCut,
    Truth_TauMatch_reldPtvsLandEtaCut,
    Truth_TauMatch_dPtvsPhiandEtaCut,
    Truth_TauMatch_dEvsPhiandEtaCut,
    Truth_TauMatch_dphivsPhiandEtaCut,
    Truth_TauMatch_dthetavsPhiandEtaCut,
    Truth_TauMatch_reldPtvsPhiandEtaCut,
    Truth_PionDr,
    Truth_PionDrAll;


  std::map<unsigned int,unsigned int> JAKIDtoIndex;
 

};

#endif
