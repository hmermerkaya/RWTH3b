// -*- C++ -*-
//
// Package:    ControlSample_SkimCuts
// Class:      ControlSample_SkimCuts
// 
/**\class ControlSample_SkimCuts ControlSample_SkimCuts.cc SkimmingTools/SkimmingCuts/src/ControlSample_SkimCuts.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian Nugent
//         Created:  Thu March 30 10:52:30 CET 2012
//
//


// system include files
#include <memory>
#include "TLorentzVector.h"
#include "TVector3.h"

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
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
//
// class declaration
//

class ControlSample_SkimCuts : public edm::EDFilter {
   public:
      explicit ControlSample_SkimCuts(const edm::ParameterSet&);
      ~ControlSample_SkimCuts();

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
  edm::InputTag muonsTag_;
  edm::InputTag pfjetsTag_;


  double MuonPtCut_;
  double MuonEtaCut_;
  double JetEtaCuts_;
  double JetPtCuts_;
  double dphiMuJet_;



  int cnt_;
  int cntFound_;

};


ControlSample_SkimCuts::ControlSample_SkimCuts(const edm::ParameterSet& iConfig):
  muonsTag_(iConfig.getParameter<edm::InputTag>( "muonsTag" )),
  pfjetsTag_( iConfig.getParameter<edm::InputTag>( "pfjetsTag" ) ),
  MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
  MuonEtaCut_(iConfig.getParameter<double>("MuonEtaCut")),
  JetEtaCuts_(iConfig.getParameter<double>("JetEtaCut")),
  JetPtCuts_(iConfig.getParameter<double>("JetPtCut")),
  dphiMuJet_(iConfig.getParameter<double>("dphiMuJet"))
{
   //now do what ever initialization is needed

}


ControlSample_SkimCuts::~ControlSample_SkimCuts()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool ControlSample_SkimCuts::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  iEvent_=&iEvent;
  cnt_++;
  
  //// Muons
  std::vector<TLorentzVector> Muon;
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel(muonsTag_,  muonCollection);
  unsigned int NMuonsPass(0), Muon_index(0);
  for(reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon!= muonCollection->end(); ++iMuon, Muon_index++){
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    if(RefMuon->p4().Pt() > MuonPtCut_){
      if(fabs(RefMuon->p4().Eta())<MuonEtaCut_){
	if(RefMuon->isGlobalMuon() && RefMuon->isStandAloneMuon()){
	  if(RefMuon->globalTrack()->normalizedChi2()<10.0){
	    if(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits()>0){
	      const reco::MuonIsolation  Iso05 = RefMuon->isolationR05();
	      if(((Iso05.emEt + Iso05.hadEt + Iso05.sumPt)/RefMuon->p4().Pt())<0.5){
		NMuonsPass++;
		Muon.push_back(TLorentzVector(RefMuon->p4().Px(),RefMuon->p4().Py(),RefMuon->p4().Pz(),RefMuon->p4().E()));
	      }
	    }
	  }
	}
      }
    }
  }
  if(NMuonsPass==0) return false;
  if(NMuonsPass>=2){cntFound_++; return true;}
 
  // try get matching
  edm::Handle<reco::PFJetCollection> JetCollection;
  iEvent.getByLabel(pfjetsTag_,  JetCollection);
  for(reco::PFJetCollection::size_type iPFJet = 0; iPFJet < JetCollection->size(); iPFJet++) {
    reco::PFJetRef PFJet(JetCollection, iPFJet);
    if(fabs(PFJet->p4().Eta())<JetEtaCuts_ && JetPtCuts_<PFJet->p4().Pt()){
      for(unsigned int i=0; i<Muon.size();i++){
	TVector3 Jet(PFJet->p4().Px(),PFJet->p4().Py(),PFJet->p4().Pz());
	std::cout << "Jet Loop angle " << fabs(Jet.DeltaPhi(Muon.at(i).Vect())) << std::endl;
	if(fabs(Jet.DeltaPhi(Muon.at(i).Vect()))>dphiMuJet_){cntFound_++; return true;}
      }
    }
  }
  return false;
}




// ------------ method called once each job just before starting event loop  ------------
void 
ControlSample_SkimCuts::beginJob()
{
  cnt_ =0;
  cntFound_ =0;

 }

// ------------ method called once each job just after ending the event loop  ------------
void 
ControlSample_SkimCuts::endJob() {
  std::cout << "--> [ControlSample_SkimCuts] found at least "<< cntFound_ <<" candidate(s). Efficiency: " <<  cntFound_ << "/" << cnt_ << " = " << ((float)cntFound_)/((float)cnt_) << "%"<< std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
ControlSample_SkimCuts::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ControlSample_SkimCuts::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ControlSample_SkimCuts::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ControlSample_SkimCuts::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ControlSample_SkimCuts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ControlSample_SkimCuts);
