#include "TriggerFilter/Filter/interface/TriggerFilterInfoProducer.h"

TriggerFilterInfoProducer::TriggerFilterInfoProducer( const edm::ParameterSet& pset ) :
  MyTriggerHelper()
{
  produces<std::vector<std::string> >("TriggerFilterInfoList").setBranchAlias("TriggerFilterInfoList");
}

TriggerFilterInfoProducer::~TriggerFilterInfoProducer(){
}

void TriggerFilterInfoProducer::beginJob(){
  return;
}

void TriggerFilterInfoProducer::produce( edm::Event& e, const edm::EventSetup& iSetup){
  std::auto_ptr<std::vector<std::string> > TriggerFilterInfoList(new std::vector<std::string>());
  *TriggerFilterInfoList = MyTriggerHelper.GetTriggerList();    
  e.put(TriggerFilterInfoList,"TriggerFilterInfoList");  
  return ;
}  

void TriggerFilterInfoProducer::endRun( const edm::Run& r, const edm::EventSetup& iSetup){
  return;
}


void TriggerFilterInfoProducer::endJob(){
  return;
}

