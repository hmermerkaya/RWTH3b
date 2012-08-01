#include "RecoTauTag/KinematicTau/interface/KinematicTauSkim.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

KinematicTauSkim::KinematicTauSkim(const edm::ParameterSet& iConfig):
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  KinematicFitTauTag_(iConfig.getParameter<edm::InputTag>("KinematicFitTauTag")),
  minTau_( iConfig.getUntrackedParameter<unsigned int>( "minTau", 1 ) )//filter returns true if more than minTau_ taus are found
{
}

KinematicTauSkim::~KinematicTauSkim(){
}

bool KinematicTauSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool filterValue = false;
  cnt_++;

  edm::Handle<SelectedKinematicDecayCollection> KinematicFitTaus;
  iEvent.getByLabel(KinematicFitTauTag_,KinematicFitTaus);
  
  unsigned int ntaus(0);  
  for(SelectedKinematicDecayCollection::const_iterator kinFitTau=KinematicFitTaus->begin();kinFitTau!=KinematicFitTaus->end();kinFitTau++){
    bool passed=true;
    for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
      if(kinFitTau->discriminators().count(*discr)>0){
	if(!(kinFitTau->discriminators().find(*discr)->second))passed=false;
      }
    }
    if(passed)ntaus++;		
  }
  edm::LogInfo("KinematicTauSkim")<<"ntaus " << ntaus << " required " << minTau_ << " inputsize " << KinematicFitTaus->size();
  if(ntaus >= minTau_){
    cntFound_++;//found at least minTau_ refitted tau(s)
    filterValue = true;
  }
  
  return filterValue;
}

void KinematicTauSkim::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;
}

void KinematicTauSkim::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("KinematicTauSkim")<<"--> [KinematicTauSkim] asks for >= "<<minTau_<<" tau(s) per event passing the provided discriminators. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauSkim);
