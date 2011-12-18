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

    std::string TrigerTemplateDATA_1 = "HLT_IsoMu12_LooseIsoPFTau10";  // reco-1
    std::string TrigerTemplateDATA_2 = "HLT_IsoMu15_LooseIsoPFTau15";  // reco-1

    std::string TrigerTemplateDATA_6 = "HLT_IsoMu15_LooseIsoPFTau15";  // reco-6
    std::string TrigerTemplateMC = "HLT_IsoMu15_LooseIsoPFTau15";




    std::string temp1 = "HLT_IsoMu";
    std::string temp2 = "_LooseIsoPFTau";

    int ntrigs=TriggerResults_->size();
    for (int itrig = 0; itrig < ntrigs; ++itrig) {
      //std::cout << trigNames.triggerName(itrig) << std::endl;
      
      
      std::string str1;
      std::string str2;

      str1 = trigNames.triggerName(itrig).substr(0,9);
      //   std::cout<<"size " <<trigNames.triggerName(itrig).size() <<std::endl;
      if(str1 == temp1 && trigNames.triggerName(itrig).size() > 20){

      str2 = trigNames.triggerName(itrig).substr(11,14);

      //   std::cout<<str1 <<"  "<<str2 <<std::endl;

      if(str2 == temp2){
	std::cout<<"Coincide: " <<str1 <<"   "<<str2 <<"  " <<trigNames.triggerName(itrig) <<std::endl;
 	accept = (*TriggerResults_).accept(itrig);
 	passedTriggerName = trigNames.triggerName(itrig);
	

//       if(trigNames.triggerName(itrig)  == "HLT_IsoMu15_LooseIsoPFTau15_v8"){ 
// 	// HLT_IsoMu15_LooseIsoPFTau15_v8
// 	// HLT_IsoMu15_LooseIsoPFTau20_v6
// 	// HLT_IsoMu15_TightIsoPFTau20_v6
// 	accept = (*TriggerResults_).accept(itrig);
// 	passedTriggerName = trigNames.triggerName(itrig);
//       }
      }
      }
    }
    
  }else{
    std::cout<< "WARNING! Trigger results is not valid! "<<std::endl;
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
