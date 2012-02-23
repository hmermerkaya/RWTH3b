// -*- C++ -*-
//
// Package:    Filter
// Class:      Filter
// 
/**\class Filter Filter.cc TriggerFilter/Filter/src/Filter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Cherepanov
//         Created:  Wed Nov 23 16:24:12 CET 2011
// $Id: Filter.cc,v 1.1 2011/12/18 15:12:01 cherepan Exp $
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
#include "DataFormats/Common/interface/TriggerResults.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class declaration
//

class Filter : public edm::EDFilter {
   public:
      explicit Filter(const edm::ParameterSet&);
      ~Filter();

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
  edm::Event * iEvent_;
  std::string processName_; // process name of (HLT) process for which to get HLT configuration
  HLTConfigProvider hltConfig_;/// The instance of the HLTConfigProvider as a data member
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
Filter::Filter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


Filter::~Filter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool accept = false;
  bool Muon_Tau = false;
  bool Muon_Eta_Tau = false;
  bool Muon_Tau_accept = false;
  bool Muon_Eta_Tau_accept = false;

  iEvent_ = &iEvent;
  cnt_++;
  
  //   char TrigerTemplate[28];
  //   sprintf(TrigerTemplate,"HLT_IsoMu%d_LooseIsoPFTau%d",15,15);
  
  
  
  edm::InputTag trigResultsTag("TriggerResults","","HLT");
  
  
  edm::Handle<edm::TriggerResults > TriggerResults_;
  iEvent_->getByLabel(trigResultsTag, TriggerResults_);
  
  
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*TriggerResults_);   
  std::string passedTriggerName;
  if (TriggerResults_.isValid()) {
    std::string checkMuon = "HLT_IsoMu";
    std::string checkTau = "_LooseIsoPFTau";
    std::string checkEta = "eta2p1";




    //-----------------  Triggers to be checked ----------------------    
    //HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1
    //HLT_IsoMu15_LooseIsoPFTau15_v9
    //---------------------------------------------------------------
    //if both triggers found the result of
    //HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1  is taken
    //
    //----------------------------------------------------------------
    int ntrigs=TriggerResults_->size();
    for (int itrig = 0; itrig < ntrigs; ++itrig) {
      //std::cout << trigNames.triggerName(itrig) << std::endl;
      
      std::string Muon;
      std::string Tau;
      std::string Eta;
      
      Muon = trigNames.triggerName(itrig).substr(0,9);
      // std::cout<<"size " <<trigNames.triggerName(itrig).size() <<std::endl;
      if(Muon == checkMuon && trigNames.triggerName(itrig).size() > 20){
	
	Tau = trigNames.triggerName(itrig).substr(11,14);
	Eta = trigNames.triggerName(itrig).substr(12,6);
	//std::cout<<Muon <<"  "<<Tau <<"  "<<Eta<<std::endl;
	if(Tau == checkTau){
	  //std::cout<<"HLT_IsoMu15_LooseIsoPFTau15  FOUND!  " <<Muon <<"   "<<Tau <<"  " <<trigNames.triggerName(itrig) <<std::endl;
	  Muon_Tau = true;
	  Muon_Tau_accept = (*TriggerResults_).accept(itrig);
	  passedTriggerName = trigNames.triggerName(itrig);
	  
	}else if(Eta == checkEta){
	  //std::cout<<"HLT_IsoMu15_LooseIsoPFTau15  Not  FOUND!  " <<std::endl;
	  Tau =  trigNames.triggerName(itrig).substr(18,14);
	  if(Tau == checkTau){
	    Muon_Eta_Tau = true;
	    Muon_Eta_Tau_accept = (*TriggerResults_).accept(itrig);
	    passedTriggerName = trigNames.triggerName(itrig);
	  }
	}
      }
    }
    
  }else{
    std::cout<< "WARNING! Trigger results is not valid! "<<std::endl;
  }
  if(Muon_Tau && !Muon_Eta_Tau){
    accept = Muon_Tau_accept;
    // std::cout<<"Only first trigger found; Return its result  " << Muon_Tau_accept <<std::endl;

  }
  if(!Muon_Tau && Muon_Eta_Tau){

    accept = Muon_Eta_Tau_accept;
    //  std::cout<<"Only second trigger found; Return its result  " << Muon_Eta_Tau_accept <<std::endl;

  }
  if(Muon_Tau && Muon_Eta_Tau){

    accept = Muon_Eta_Tau_accept;
    //  std::cout<<"Both triggers found; Return result of the second  " << Muon_Tau_accept << " <>  "<< Muon_Eta_Tau_accept <<std::endl;

  }

  if(!Muon_Tau && !Muon_Eta_Tau){

    accept = false;
    //std::cout<<"No triggers found; Return FALSE" <<std::endl;

  }

  if(accept) cntFound_++;
  
  return accept;
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
Filter::beginJob()
{
  cnt_=0;
  cntFound_=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Filter::endJob() {
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  printf("-> [Filter] Trigger efficiency %d/%d = %.2f%%\n", cntFound_, cnt_, ratio*100.0);


}


// ------------ method called when starting to processes a run  ------------
bool 
Filter::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup){ 
  std::string processName = "HLT";
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName,changed)) {
    // if init returns TRUE, initialisation has succeeded!
     printf(" %d \n",changed);
    if (changed) {
      const std::string &  DSName = hltConfig_.datasetName(1);
      std::cout << DSName << std::endl;
      const std::vector< std::vector< std::string > > & AllDSName=hltConfig_.datasetContents(); 
      const std::vector< std::string > & TriggNames=hltConfig_.triggerNames();


      for(std::vector< std::string >::const_iterator iName = TriggNames.begin(); iName !=TriggNames.end(); ++iName ){

			std::cout << (*iName) << std::endl;
      }
     // The HLT config has actually changed wrt the previous Run, hence rebook your
     // histograms or do anything else dependent on the revised HLT config
      std::cout<<"HLTSelector:: tableName of process "<<processName<<" is "<<hltConfig_.tableName()<<std::endl;
//  printf("succes \n");
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    std::cout<<"MyAnalyzer" << " HLT config extraction failure with process name " << processName_<<std::endl;
    // In this case, all access methods will return empty values!
  }

  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
Filter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
Filter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
Filter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Filter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Filter);
