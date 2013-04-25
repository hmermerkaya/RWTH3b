#include "RecoTauTag/KinematicTau/interface/Tau_JAKID_Filter.h"

#include "Validation/EventGenerator/interface/TauDecay.h"
#include "RecoTauTag/KinematicTau/interface/TauDecay_CMSSWReco.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TLorentzVector.h"

Tau_JAKID_Filter::Tau_JAKID_Filter(const edm::ParameterSet& iConfig):
  JAKID_( iConfig.getParameter< std::vector<int> >("jakid") ),
  nprongs_( iConfig.getParameter< std::vector<int> >("nprongs") ),
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
  unsigned int ntaus=0;  
  unsigned int n3prong=0;
  if(genParticles.isValid()){
    for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
     const reco::GenParticle mytau=(*itr);
     if(fabs(mytau.pdgId())==fabs(PdtPdgMini::tau_minus)){
       TauDecay_CMSSWReco TD;
       unsigned int jak_id, TauBitMask;
       TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
       if(isTruthTauInAcceptance(mytau,TauPtMin_,TauEtaMax_)){
	 for(unsigned int i=0;i<JAKID_.size() && i<nprongs_.size();i++){
	   //std::cout << "JAKID " << JAKID_.at(i) << " " << jak_id << " " << nprongs_.at(i) << " " << TD.nProng(TauBitMask) << std::endl;
	   if(JAKID_.at(i)==(int)jak_id && (int)TD.nProng(TauBitMask)==nprongs_.at(i)){ n3prong++;}
	   if(JAKID_.at(i)==(int)jak_id){ntaus++;}
	 }
       }
     }
    }
  }
  if(ntaus==1 && n3prong==1)filterValue=true;  
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
  std::cout << "Tau_JAKID_Filter" <<"--> [Tau_JAKID_Filter]  Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%" << std::endl;
  for(unsigned int i=0; i<JAKID_.size() && i<nprongs_.size();i++){
    std::cout << "JAKID " << JAKID_.at(i) << " " << nprongs_.at(i) << std::endl;
  }
}

// This function is duplicated... needs to be fixed!!
bool Tau_JAKID_Filter::isTruthTauInAcceptance(const reco::GenParticle &cand,double &TauPtMin,double &TauEtaMax){
  if(fabs(cand.pdgId())!=fabs(PdtPdgMini::tau_minus)) return false;
  //if(cand.status()!=2) return false; // require tau that is: 2) a particle after parton showering and ISR/FSR
  for (unsigned int i=0; i< cand.numberOfDaughters(); i++){
    if(cand.daughter(i)->pdgId()==cand.pdgId()) return false;

  }

  TLorentzVector mc_tau(cand.p4().Px(),cand.p4().Py(),cand.p4().Pz(),cand.p4().E());
  //if(tau.Pt()>TauPtMin_ && fabs(tau.Eta())<TauEtaMax_)return true; // require tau within Pt and |eta| acceptance
  TauDecay_CMSSWReco TD;
  unsigned int jak_id, TauBitMask;
  TD.AnalyzeTau(&cand,jak_id,TauBitMask);
  std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
  TLorentzVector mc_nu(0,0,0,0);
  for(unsigned int j=0;j<DecayProd.size();j++){
    if(fabs(DecayProd.at(j)->pdgId())==PdtPdgMini::nu_tau){
      mc_nu.SetPxPyPzE(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
    }
  }
  TLorentzVector vis=mc_tau-mc_nu;
  if(vis.Pt()>TauPtMin && fabs(vis.Eta())<TauEtaMax)return true; // require tau within Pt and |eta| acceptance                 
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Tau_JAKID_Filter);
