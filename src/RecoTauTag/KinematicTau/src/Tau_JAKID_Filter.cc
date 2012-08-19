#include "RecoTauTag/KinematicTau/interface/Tau_JAKID_Filter.h"

#include "Validation/EventGenerator/interface/TauDecay.h"
#include "RecoTauTag/KinematicTau/interface/TauDecay_CMSSWReco.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TLorentzVector.h"

Tau_JAKID_Filter::Tau_JAKID_Filter(const edm::ParameterSet& iConfig):
  JAKID_( iConfig.getParameter< std::vector<int> >("jakid") ),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  TauPtMin_( iConfig.getParameter<double>("TauPtMin")),
  TauEtaMax_( iConfig.getParameter<double>("TauEtaMax"))
{
}

Tau_JAKID_Filter::~Tau_JAKID_Filter(){
}

bool Tau_JAKID_Filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool filterValue = false;
  cnt_++;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(gensrc_, genParticles);
  
  if(genParticles.isValid()){
    for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
      const reco::GenParticle mytau=(*itr);
      if(isTruthTauInAcceptance(mytau)){
	TauDecay_CMSSWReco TD;
	unsigned int jak_id, TauBitMask;
	TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	for(unsigned int i=0;i<JAKID_.size();i++){
	  if(((unsigned int)JAKID_.at(i))==jak_id){ filterValue=true; break;}
	}
	if(filterValue) break;
      }
    }
  }
  
  if(filterValue){
    cntFound_++;
  }
  
  return filterValue;
}

void Tau_JAKID_Filter::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;
}

void Tau_JAKID_Filter::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("Tau_JAKID_Filter")<<"--> [Tau_JAKID_Filter]  Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

// This function is duplicated... needs to be fixed!!
bool Tau_JAKID_Filter::isTruthTauInAcceptance(const reco::GenParticle &cand){
  if(fabs(cand.pdgId())!=fabs(PdtPdgMini::tau_minus)) return false;
  if(cand.status()!=2) return false; // require tau that is: 2) a particle after parton showering and ISR/FSR
  TLorentzVector tau(cand.p4().Px(),cand.p4().Py(),cand.p4().Pz(),cand.p4().E());
  if(tau.Pt()>TauPtMin_ && fabs(tau.Eta())<TauEtaMax_)return true; // require tau within Pt and |eta| acceptance
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Tau_JAKID_Filter);
