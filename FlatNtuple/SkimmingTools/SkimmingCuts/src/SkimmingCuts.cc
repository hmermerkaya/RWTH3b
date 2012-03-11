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
// $Id$
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
  bool PFtausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  edm::Event * iEvent_;
  edm::InputTag hpsTauProducer_;


  double MuonPtCut_;
  bool MuonIsGlo_;
  double NMuons_;
  double MuonEtaCut_;
  double PFTauPtCut_;
  double PFTauEtaCut_;


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
PFTauEtaCut_( iConfig.getParameter<double>("PFTauEtaCut") )
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
  bool AcceptPFTau = PFtausCuts(iEvent, iSetup);

  if(AcceptMuon && AcceptPFTau){
    pass = true;
    cntFound_++;
  }
   return pass;
}

bool 
SkimmingCuts::MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;
  bool nMuonsCut = false;
  bool MuonPtCut = false;
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel("muons",  muonCollection);

  int Muon_index =0;
  reco::MuonRef HighestPtMuonRef;
  double Pt = 0;
  for(reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon!= muonCollection->end(); ++iMuon, Muon_index++){
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    if(RefMuon->p4().Pt() > Pt){
      Pt = RefMuon->p4().Pt();
      HighestPtMuonRef = RefMuon;
    }
  }
  if(HighestPtMuonRef.isNonnull()){
    nMuon_++;
    if(HighestPtMuonRef->p4().Pt() > MuonPtCut_){
      if(fabs(HighestPtMuonRef->p4().Eta()) < MuonEtaCut_){
	if(HighestPtMuonRef->isGlobalMuon() == MuonIsGlo_){
	  MuonPtCut = true;
	  nMuonPass_++;
	  //  std::cout<<"Muon:  "<<"Pt: " <<HighestPtMuonRef->p4().Pt() << "   global: " << HighestPtMuonRef->isGlobalMuon()  << "   eta: " << HighestPtMuonRef->p4().Eta()<<std::endl;
	}
      }
    }
  }
   if(nMuonPass_ == NMuons_){
     nMuonsCut = true;
   }

  pass = nMuonsCut*nMuonsCut;

  return pass;
}

bool 
SkimmingCuts::PFtausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool pass = false;
  edm::Handle<std::vector<reco::PFTau> > HPStaus;
  iEvent_->getByLabel(hpsTauProducer_, HPStaus);
  
  reco::PFTauRef HighestPtPFTauCandidate;
  double Pt = 0;

  for ( unsigned iPFTau = 0; iPFTau < HPStaus->size(); ++iPFTau ) {
    
    reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);

    if(HPStauCandidate->p4().Pt() > Pt){
      Pt = HPStauCandidate->p4().Pt();
      HighestPtPFTauCandidate = HPStauCandidate;
    }
  }

  if(HighestPtPFTauCandidate.isNonnull()){
    nPFTaus_++;
    if(HighestPtPFTauCandidate->p4().Pt() > PFTauPtCut_){
      if(fabs(HighestPtPFTauCandidate->p4().Eta()) < PFTauEtaCut_){
	  pass = true;
	  nPFTausPass_++;
	  // std::cout<<"PFTau:  "<<"Pt: " <<HighestPtPFTauCandidate->p4().Pt() <<"   eta: " << HighestPtPFTauCandidate->p4().Eta()<<std::endl;

      }
    }
  }


  return pass;
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


   std::cout<<"Starting preselection ...  " <<std::endl;
 }

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimmingCuts::endJob() {

  float ratioMuon = 0.0;
  if(nMuon_!=0) ratioMuon=(float)nMuonPass_/nMuon_;

  float ratioPFTau = 0.0;
  if(nPFTaus_!=0) ratioPFTau=(float)nPFTausPass_/nPFTaus_;


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
