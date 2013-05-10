#include "TriggerFilter/Filter/interface/TriggerFilter.h"


TriggerFilter::TriggerFilter(const edm::ParameterSet& iConfig) :
  MyTriggerHelper(),
  doTauplusXTrigger_(iConfig.getUntrackedParameter("doTauplusXTrigger",(bool)(true))), // set as default to keep consistent with orig. version
  doMuonTrigger_(iConfig.getUntrackedParameter("doMuonTrigger",(bool)(false))),
  doElectronTrigger_(iConfig.getUntrackedParameter("doElectronTrigger",(bool)(false))),
  doTauplusMETTrigger_(iConfig.getUntrackedParameter("doTauplusMETTrigger",(bool)(false)))
{
   //now do what ever initialization is needed

}


TriggerFilter::~TriggerFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  MyTriggerHelper.Reset();
  bool TauplusXTrigger(false),MuonTrigger(false),ElectronTrigger(false),TauplusMETTrigger(false);
  if(doTauplusXTrigger_)   TauplusXTrigger= FilteronTauplusXTrigger(iEvent);
  if(doMuonTrigger_)       MuonTrigger=FilteronMuonTrigger(iEvent);
  if(doElectronTrigger_)   ElectronTrigger=FilteronElectronTrigger(iEvent);
  if(doTauplusMETTrigger_) TauplusMETTrigger=FilteronTauplusMETTrigger(iEvent);
  //  if(TauplusXTrigger || MuonTrigger || ElectronTrigger || TauplusMETTrigger ) return true;
  return true;
  //return false;
}

bool TriggerFilter::FilteronTauplusXTrigger(edm::Event& iEvent){
  // Vladimir Cherepanov original TriggerFilter
  bool accept = false;
  bool Muon_Tau = false;
  bool Muon_Eta_Tau = false;
  bool Muon_Tau_accept = false;
  bool Muon_Eta_Tau_accept = false;

  iEvent_ = &iEvent;
  cnt_++;
  
  
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
      
      std::string TriggerNa;

//       Muon = trigNames.triggerName(itrig).substr(0,9);
//       TriggerNa=trigNames.triggerName(itrig).substr(0,34);
//       // std::cout<<"size " <<trigNames.triggerName(itrig).size() <<std::endl;
//       if(Muon == checkMuon && trigNames.triggerName(itrig).size() > 20){
	
// 	Tau = trigNames.triggerName(itrig).substr(11,14);
// 	Eta = trigNames.triggerName(itrig).substr(12,6);
// 	//	std::cout<<Muon <<"  "<<Tau <<"  "<<Eta<<"  all  " << Muon + Eta + Tau<<std::endl;
// 	//	std::cout<<"TriggerNa  "<< TriggerNa<<std::endl;
// 	if(Tau == checkTau){
// 	  //std::cout<<"HLT_IsoMu15_LooseIsoPFTau15  FOUND!  " <<Muon <<"   "<<Tau <<"  " <<trigNames.triggerName(itrig) <<std::endl;
// 	  Muon_Tau = true;
// 	  Muon_Tau_accept = (*TriggerResults_).accept(itrig);
// 	  passedTriggerName = trigNames.triggerName(itrig);
	  
// 	}else if(Eta == checkEta){
// 	  //std::cout<<"HLT_IsoMu15_LooseIsoPFTau15  Not  FOUND!  " <<std::endl;
// 	  Tau =  trigNames.triggerName(itrig).substr(18,14);
// 	  if(Tau == checkTau){
// 	    Muon_Eta_Tau = true;
// 	    Muon_Eta_Tau_accept = (*TriggerResults_).accept(itrig);
// 	    passedTriggerName = trigNames.triggerName(itrig);
// 	  }
// 	}
//       }


      if(trigNames.triggerName(itrig).substr(0,34) == "HLT_IsoMu18_eta2p1_LooseIsoPFTau20" or trigNames.triggerName(itrig).substr(0,34) == "HLT_IsoMu17_eta2p1_LooseIsoPFTau20" ){

	passedTriggerName = trigNames.triggerName(itrig);
	accept = (*TriggerResults_).accept(itrig);
      }

    }
    
  }else{
    std::cout<< "WARNING! Trigger results is not valid! "<<std::endl;
  }
//   if(Muon_Tau && !Muon_Eta_Tau){
//     accept = Muon_Tau_accept;
//     //   std::cout<<"Only first trigger found; Return its result  " << Muon_Tau_accept <<std::endl;

//   }
//   if(!Muon_Tau && Muon_Eta_Tau){

//     accept = Muon_Eta_Tau_accept;
//     //  std::cout<<"Only second trigger found; Return its result  " << Muon_Eta_Tau_accept <<std::endl;

//   }
//   if(Muon_Tau && Muon_Eta_Tau){

//     accept = Muon_Eta_Tau_accept;
//     //   std::cout<<"Both triggers found; Return result of the second  " << Muon_Tau_accept << " <>  "<< Muon_Eta_Tau_accept <<std::endl;

//   }

//   if(!Muon_Tau && !Muon_Eta_Tau){

//     //  accept = false;
//     //std::cout<<"No triggers found; Return FALSE" <<std::endl;

//   }



  if(accept) cntFound_++;
  // Store Trigger name for Trigger Helper
  MyTriggerHelper.AddTrigger(passedTriggerName);
  //  std::cout<<"                  passedTriggerName                     "<<passedTriggerName<<std::endl;
  return accept;
}


bool TriggerFilter::FilteronMuonTrigger(edm::Event& iEvent){

  bool accept = false;
  bool Muon = false;

  iEvent_ = &iEvent;
  cnt_++;
  
  
  edm::InputTag trigResultsTag("TriggerResults","","HLT");
  edm::Handle<edm::TriggerResults > TriggerResults_;
  iEvent_->getByLabel(trigResultsTag, TriggerResults_);
  
  
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*TriggerResults_);   
  std::string passedTriggerName;
  if (TriggerResults_.isValid()) {
    std::string checkMuon = "HLT_IsoMu24_v";
    int ntrigs=TriggerResults_->size();
    for (int itrig = 0; itrig < ntrigs; ++itrig) {
      //std::cout << trigNames.triggerName(itrig) << std::endl;
      std::string Muon;
      Muon = trigNames.triggerName(itrig).substr(0,13);
      if(Muon == checkMuon && trigNames.triggerName(itrig).size() < 20){
	accept = (*TriggerResults_).accept(itrig);
	//	std::cout<<"HLT_IsoMu24  FOUND!  " <<Muon <<"   " <<"  " <<trigNames.triggerName(itrig) <<std::endl;
	passedTriggerName = trigNames.triggerName(itrig);
      }
    }
  }
  
  
  if(accept) cntFound_++;
  // Store Trigger name for Trigger Helper
  MyTriggerHelper.AddTrigger(passedTriggerName);
  return accept;


}
bool TriggerFilter::FilteronElectronTrigger(edm::Event& iEvent){
  return false;
}
bool TriggerFilter::FilteronTauplusMETTrigger(edm::Event& iEvent){
  return false;
}



// ------------ method called once each job just before starting event loop  ------------
void 
TriggerFilter::beginJob()
{
  cnt_=0;
  cntFound_=0;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerFilter::endJob() {
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  printf("-> [TriggerFilter] Trigger efficiency %d/%d = %.2f%%\n", cntFound_, cnt_, ratio*100.0);


}








// ------------ method called when starting to processes a run  ------------
bool 
TriggerFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TriggerFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TriggerFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TriggerFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}







// // // ------------ method called when starting to processes a run  ------------
// // bool 
// // TriggerFilter::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup){ 
// // //   std::string processName = "HLT";
// // //   bool changed(true);
// // //   std::cout<<"-------------------- >>>>>>>>>>>>>> TriggerFilter debug 1"<<std::endl;
// // //   if (hltConfig_.init(iRun,iSetup,processName,changed)) {
// // //     // if init returns TRUE, initialisation has succeeded!
// // //     //    printf(" %d \n",changed);
// // //     if (changed) {
// // //       const std::string &  DSName = hltConfig_.datasetName(1);
// // //       //  std::cout << DSName << std::endl;
// // //       //const std::vector< std::vector< std::string > > & AllDSName=hltConfig_.datasetContents(); 
// // //       const std::vector< std::string > & TriggNames=hltConfig_.triggerNames();

      
// // //       //       for(std::vector< std::string >::const_iterator iName = TriggNames.begin(); iName !=TriggNames.end(); ++iName ){
      
// // //       // 	//		std::cout << (*iName) << std::endl;
// // //       //       }
// // //       // The HLT config has actually changed wrt the previous Run, hence rebook your
// // //       // histograms or do anything else dependent on the revised HLT config
// // //       std::cout<<"HLTSelector:: tableName of process "<<processName<<" is "<<hltConfig_.tableName()<<std::endl;
      
// // //       //  printf("succes \n");
// // //     }
// // //   } else {
// // //     // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
// // //     // with the file and/or code and needs to be investigated!
// // //     std::cout<<"MyAnalyzer" << " HLT config extraction failure with process name " << processName_<<std::endl;
// // //     // In this case, all access methods will return empty values!
// // //   }
// //   std::cout<<"-------------------- >>>>>>>>>>>>>> TriggerFilter debug 2"<<std::endl;
// //   return true;
// //   std::cout<<"-------------------- >>>>>>>>>>>>>> TriggerFilter debug 3"<<std::endl;
  
// // }
// // // ------------ method called when ending the processing of a run  ------------
// // bool 
// // TriggerFilter::endRun(edm::Run&, edm::EventSetup const&)
// // {
// // std::cout<<"-------------------- >>>>>>>>>>>>>> TriggerFilter End Run"<<std::endl;
// //   return true;
// // }

// // ------------ method called when starting to processes a luminosity block  ------------
// bool 
// TriggerFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
// {
//   return true;
// }

// // ------------ method called when ending the processing of a luminosity block  ------------
// bool 
// TriggerFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
// {
//   return true;
// }

// // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void
// TriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }


