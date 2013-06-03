// -*- C++ -*-
//
// Package:    SkimmingCuts
// Class:      SkimmingCuts
// 
/**\class SkimmingCuts SkimmingCuts.cc SkimmingTools/SkimmingCuts/src/SkimmingCuts.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Cherepanov
//         Created:  Thu Feb 23 18:53:30 CET 2012
// $Id: SkimmingCuts.cc,v 1.3 2013/05/10 09:54:25 cherepan Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/MuonReco/interface/MuonIsolation.h>
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include <TMath.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"

#include "TVectorT.h"
#include "TH1.h"
#include "TLorentzVector.h"


//
// class declaration
//

class SkimmingCuts : public edm::EDFilter {
   public:
      explicit SkimmingCuts(const edm::ParameterSet&);
      ~SkimmingCuts();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  bool MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool KFitTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool PFTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  reco::PFTauRef getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &UnmodifiedTau, unsigned int &match);
  double DeltaPhi(double phi1, double phi2);
  edm::Event * iEvent_;
  edm::InputTag hpsTauProducer_;


  double MuonPtCut_;
  bool MuonIsGlo_;
  double NMuons_;
  double MuonEtaCut_;
  double PFTauPtCut_;
  double PFTauEtaCut_;

  double ElectronPtCut_;
  double ElectronEtaCut_;
  edm::InputTag KinFitAdvanced_;

  int nMuon_;
  int nMuonPass_;
  int nPFTaus_;
  int nPFTausPass_;


  int cnt_;
  int cntFound_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SkimmingCuts::SkimmingCuts(const edm::ParameterSet& iConfig):
hpsTauProducer_( iConfig.getParameter<edm::InputTag>( "hpsTauProducer" ) ),
MuonPtCut_( iConfig.getParameter<double>("MuonPtCut") ),
MuonIsGlo_( iConfig.getParameter<bool>("MuonIsGlobal") ),
NMuons_( iConfig.getParameter<double>("NMuons") ),
MuonEtaCut_( iConfig.getParameter<double>("MuonEtaCut") ),
PFTauPtCut_( iConfig.getParameter<double>("PFTauPtCut") ),
PFTauEtaCut_( iConfig.getParameter<double>("PFTauEtaCut") ),
ElectronPtCut_( iConfig.getParameter<double>("ElectronPtCut") ),
ElectronEtaCut_( iConfig.getParameter<double>("ElectronEtaCut") ),
KinFitAdvanced_( iConfig.getParameter<edm::InputTag>( "kinematicTausAdvanced" ) )
{
   //now do what ever initialization is needed

}


SkimmingCuts::~SkimmingCuts()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SkimmingCuts::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cnt_++;
  bool pass = false;
  iEvent_=&iEvent;
  bool AcceptMuon = MuonsCuts(iEvent, iSetup);
  bool AcceptKFTau = KFitTausCuts(iEvent, iSetup);
  bool AcceptPFTau = PFTausCuts(iEvent, iSetup);
  bool AcceptElectron = ElectronCuts(iEvent, iSetup);

    //----------- This Blos is for private production to analyse muon tau decay with KF
  if(AcceptMuon && (AcceptKFTau || AcceptElectron || AcceptPFTau)){
      pass = true;
      cntFound_++;
    }
    //----------- This Blos is for private production to analyse muon tau decay with KF






    //----------- This Blos is for private production to analyse /rho and a1 tau decays w/o KF
//     if(AcceptMuon && AcceptPFTau){
//      pass = true;
//      cntFound_++;
//    }
    //----------- This Blos is for private production to analyse /rho and a1 tau decays w/o KF


    //    std::cout<<"AcceptMuon=  "<<AcceptMuon<<" AcceptPFTau=  "<<AcceptPFTau<<" AcceptElectron= "<<AcceptElectron<<" pass =  "<<pass<<std::endl;
//   if(AcceptMuon){
//     cntFound_++;
//     pass = true;
//   } 
   return pass;
}

bool 
SkimmingCuts::MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel("muons",  muonCollection);

  int Muon_index =0;
  reco::MuonRef HighestPtMuonRef;
  double Pt = 0;
  for(reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon!= muonCollection->end(); iMuon++, Muon_index++){
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    if(RefMuon->p4().Pt() > Pt){
      Pt = RefMuon->p4().Pt();
      HighestPtMuonRef = RefMuon;
    }
  }
  if(HighestPtMuonRef.isNonnull()){
    if(HighestPtMuonRef->p4().Pt() > MuonPtCut_){
      if(fabs(HighestPtMuonRef->p4().Eta()) < MuonEtaCut_){
	if(HighestPtMuonRef->isGlobalMuon() == MuonIsGlo_){
	  pass = true;
	  //  std::cout<<"Muon:  "<<"Pt: " <<HighestPtMuonRef->p4().Pt() << "   global: " << HighestPtMuonRef->isGlobalMuon()  << "   eta: " << HighestPtMuonRef->p4().Eta()<<std::endl;
	}
      }
    }
  }

  return pass;
}

bool SkimmingCuts::ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;
  bool ElectronCut = false;
  edm::Handle< reco::GsfElectronCollection > electronCollection;
  iEvent_->getByLabel("gsfElectrons", electronCollection);
  
  int Electron_index = 0;
  reco::GsfElectronRef HighestPtElectronRef;
  double Pt = 0;
  for(reco::GsfElectronCollection::const_iterator iElectron=electronCollection->begin(); iElectron!=electronCollection->end(); iElectron++, Electron_index++){
    reco::GsfElectronRef RefElectron(electronCollection, Electron_index);
    if(RefElectron->p4().Pt()>Pt){
      Pt=RefElectron->p4().Pt();
      HighestPtElectronRef = RefElectron;
    }
  }
  if(HighestPtElectronRef.isNonnull()){
    if(HighestPtElectronRef->p4().Pt()>ElectronPtCut_){
      
      if(fabs(HighestPtElectronRef->p4().Eta())<ElectronEtaCut_){
	ElectronCut = true;
      }
    }
  }
  
  pass = ElectronCut;
  
  return pass;
}


bool 
SkimmingCuts::KFitTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;

  //========= HPS taus for matching issues                                                                                                                                                                                                  
  edm::Handle<std::vector<reco::PFTau> > HPStaus;
  iEvent.getByLabel("hpsPFTauProducer", HPStaus);
  //========= HPS taus for matching issues       
  edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection", HPSAgainstElectronsTight);
  
  edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection", HPSAgainstMuonTight);

   edm::Handle<reco::PFTauDiscriminator> HPSLooseIsoDiscrDBSumPtCorr;
   // iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr", HPSLooseIsoDiscrDBSumPtCorr);
   iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationMVA", HPSLooseIsoDiscrDBSumPtCorr);


  //================== KinematicFit Info ===================
  edm::Handle<SelectedKinematicDecayCollection> selected;
  iEvent.getByLabel(KinFitAdvanced_, selected);
  //  std::cout<<"SkimmingCuts: NKFtaus  "<<selected->size()<<std::endl; 
  //start loop over KinFit decays 
  double Pt = 0;
  unsigned int tauindex=0;
  SelectedKinematicDecay LeadingKFTau;
  bool hasLeadingKFTau=false;
  if(selected->size()!=0){
    for(SelectedKinematicDecayCollection::const_iterator decay = selected->begin(); decay != selected->end(); ++decay, tauindex++){
      SelectedKinematicDecay KFTau=(*decay);
      if(KFTau.Initial_a1_p4().Pt() > Pt){
	Pt = KFTau.Initial_a1_p4().Pt();
	LeadingKFTau = (*decay);
	hasLeadingKFTau=true;
      }
    }

    std::vector<float> TauVisible;
    TauVisible.push_back(LeadingKFTau.Initial_a1_p4().E());
    TauVisible.push_back(LeadingKFTau.Initial_a1_p4().Px());
    TauVisible.push_back(LeadingKFTau.Initial_a1_p4().Py());
    TauVisible.push_back(LeadingKFTau.Initial_a1_p4().Pz());

    unsigned int idx =0;
    reco::PFTauRef MatchedHPSTau = getMatchedHPSTau(HPStaus,TauVisible,idx);
    if(hasLeadingKFTau){
      if(LeadingKFTau.Initial_a1_p4().Pt() > PFTauPtCut_){
	if(fabs(LeadingKFTau.Initial_a1_p4().Eta()) < PFTauEtaCut_){
	  if((*HPSAgainstElectronsTight)[MatchedHPSTau]){
	    if((*HPSAgainstMuonTight)[MatchedHPSTau]){ 
	      pass = true;
	      // std::cout<<"pass1 " <<pass<<std::endl;
	      //	      const reco::PFCandidateRefVector & 	cands =MatchedHPSTau ->signalPFChargedHadrCands(); //candidates in signal cone 
	    }
	  }
	}
      }
    }
  }
  return pass;
}

bool 
SkimmingCuts::PFTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;

  //========= HPS taus for matching issues                                                                                                                                                                                                  
  edm::Handle<std::vector<reco::PFTau> > HPStaus;
  iEvent.getByLabel("hpsPFTauProducer", HPStaus);
  //========= HPS taus for matching issues       
  edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection", HPSAgainstElectronsTight);
  
  edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection", HPSAgainstMuonTight);

  edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding", HPSByDecayModeFinding);

  edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseIsolationMVA2;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationMVA2", HPSPFTauDiscriminationByLooseIsolationMVA2);
  
  edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumIsolationMVA2;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationMVA2", HPSPFTauDiscriminationByMediumIsolationMVA2);
  
  
  edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightIsolationMVA2;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightIsolationMVA2", HPSPFTauDiscriminationByTightIsolationMVA2);



  //================== KinematicFit Info ===================
  double Pt=0;
  for ( unsigned int iPFTau = 0; iPFTau < HPStaus->size(); iPFTau++ ) {
    reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);
    if(HPStauCandidate->p4().Pt() >Pt ){
      pass=false;
      Pt = HPStauCandidate->p4().Pt();
      if(HPStauCandidate->p4().Pt() > 20){
	if((*HPSByDecayModeFinding)[HPStauCandidate]){
	  if((*HPSAgainstElectronsTight)[HPStauCandidate]){
	    if((*HPSAgainstMuonTight)[HPStauCandidate]){
	      if((*HPSPFTauDiscriminationByTightIsolationMVA2)[HPStauCandidate]){
		pass=true;
	      }
	    }
	  }
	}
      }
    }
  }
  return pass;
}


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //
 // reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &Jet, unsigned int &match)
 //
 // finds HPS tau candidate for a given KinFit tau candidate
 // the closest by deltaR HPS candidate is accepted
 reco::PFTauRef SkimmingCuts::getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &UnmodifiedTau, unsigned int &match){
   TLorentzVector TauVisible;
   TauVisible.SetE(UnmodifiedTau.at(0));
   TauVisible.SetPx(UnmodifiedTau.at(1));
   TauVisible.SetPy(UnmodifiedTau.at(2));
   TauVisible.SetPz(UnmodifiedTau.at(3));
   reco::PFTauRef MatchedHPSTau;
   double deltaR = 999;
   match=0;
    for ( unsigned int iTau = 0; iTau < HPStaus->size(); ++iTau ) {
     reco::PFTauRef HPStauCandidate(HPStaus, iTau);
     double dr=sqrt( pow(DeltaPhi(HPStauCandidate->p4().Phi(),TauVisible.Phi()),2) + pow(HPStauCandidate->p4().Eta() - TauVisible.Eta(),2));
     if(dr < deltaR){
       deltaR = dr;
       match=iTau;
       MatchedHPSTau = HPStauCandidate;
     }

   }
   return MatchedHPSTau;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// double TauNtuple::DeltaPhi(double phi1, double phi2)
//
// Calculates the the difference between two phi angles
// 
double SkimmingCuts::DeltaPhi(double phi1, double phi2){
  double dphi=fabs(phi1-phi2);
  if (dphi>TMath::Pi())   dphi=2*TMath::Pi()-dphi;
  double sign=1;
  if(phi1-phi2<0) sign=-1;
  return dphi*sign;
}



// ------------ method called once each job just before starting event loop  ------------
void 
SkimmingCuts::beginJob()
{
  nMuon_=0;
  nMuonPass_=0;
  nPFTaus_=0;
  nPFTausPass_=0;
  cnt_ =0;
  cntFound_ =0;


  //   std::cout<<"Starting preselection ...  " <<std::endl;
 }

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimmingCuts::endJob() {

//   float ratioMuon = 0.0;
//   if(nMuon_!=0) ratioMuon=(float)nMuonPass_/nMuon_;

//   float ratioPFTau = 0.0;
//   if(nPFTaus_!=0) ratioPFTau=(float)nPFTausPass_/nPFTaus_;


  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;

  //  std::cout<<"NMuons: " << nMuon_ <<"   NMuonsPass: "<< nMuonPass_ <<"   Efficiency: "<< ratioMuon*100.0 <<"%"<<std::endl;
  // std::cout<<"NPFtaus: " << nPFTaus_ <<"   NPFtausPass: "<< nPFTausPass_ <<"   Efficiency: "<< ratioPFTau*100.0 <<"%"<<std::endl;


    std::cout<<"[SkimmingCuts]-->  "<<"NEvents: " << cnt_ <<"   NEventsPass: "<< cntFound_ <<"   Efficiency: "<< ratio*100.0 <<"%"<<std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
SkimmingCuts::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SkimmingCuts::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SkimmingCuts::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SkimmingCuts::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimmingCuts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SkimmingCuts);
