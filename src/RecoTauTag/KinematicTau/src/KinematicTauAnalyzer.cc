
/**
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */
#include "RecoTauTag/KinematicTau/interface/KinematicTauAnalyzer.h"
// object includes
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "Validation/EventGenerator/interface/TauDecay.h"
#include "RecoTauTag/KinematicTau/interface/TauDecay_CMSSWReco.h"
#include "TVectorT.h"

KinematicTauAnalyzer::KinematicTauAnalyzer(const edm::ParameterSet& iConfig):
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  KinematicFitTauTag_(iConfig.getParameter<edm::InputTag>("KinematicFitTauTag")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  TauMatchingDR_( iConfig.getParameter<double>("TauMatchingDR")),
  TauPtMin_( iConfig.getParameter<double>("TauPtMin")),
  TauEtaMax_( iConfig.getParameter<double>("TauEtaMax")),
  TauM_(1.77682),     // temp solution untill particle helper made
  PionM_(0.13957018), // temp solution untill particle helper made
  NuM_(0.0)           // temp solution untill particle helper made
{
  dbe = 0;
  dbe = edm::Service<DQMStore>().operator->();
}

KinematicTauAnalyzer::~KinematicTauAnalyzer(){
}

void KinematicTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  double weight=1.0;

  edm::Handle<SelectedKinematicDecayCollection> KinematicFitTaus;
  iEvent.getByLabel(KinematicFitTauTag_,KinematicFitTaus);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(gensrc_, genParticles);

  for(SelectedKinematicDecayCollection::const_iterator kinFitTau=KinematicFitTaus->begin();kinFitTau!=KinematicFitTaus->end();kinFitTau++){
    bool passed=true;
    for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
      if(kinFitTau->discriminators().count(*discr)>0){
	if(!kinFitTau->discriminators().find(*discr)->second) passed=false;
      }
    }
    if(passed){
      SelectedKinematicDecay KFTau=(*kinFitTau);
      unsigned int ambiguity=0;
      const TLorentzVector Tau=KFTau.Tau(ambiguity);
      const TLorentzVector Nu=KFTau.Neutrino(ambiguity);
      const std::vector<TLorentzVector> Pions=KFTau.Pions(ambiguity);

      const TLorentzVector Tau_orig=KFTau.InitialTauGuess(ambiguity);
      const TLorentzVector Nu_orig=KFTau.InitalNeutrinoGuess(ambiguity);
      const std::vector<TLorentzVector> Pions_orig=KFTau.InitalPions();

      reco::Vertex Pvtx=KFTau.PrimaryVertexReFitAndRotated();
      reco::Vertex Pvtx_orig=KFTau.InitalPrimaryVertexReFitAndRotated();

      reco::Vertex Secvtx=KFTau.SecondaryVertex(ambiguity);
      reco::Vertex Secvtx_orig=KFTau.InitalSecondaryVertex();

      ///////////////////////////////
      
      const SelectedKinematicParticleCollection& Particles =KFTau.particles();
      for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = Particles.begin(); iParticle != Particles.end(); ++iParticle){
	if(iParticle->name()=="tau"){
	  TVectorT<double> intauParam;
	  TVectorT<double> parameters;
	  intauParam.ResizeTo(7);
	  parameters.ResizeTo(7);
	  intauParam=iParticle->SelectedKinematicParticle::input_parameters();
	  parameters=iParticle->SelectedKinematicParticle::parameters();
	  for(unsigned int i=0;i<7;i++){
	    std::cout << i << " Input: " << intauParam[i] << " Final " << parameters[i] << std::endl;
	  }
	}
      }
         std::cout << "Tau (" <<  Tau.Px() << ","  << Tau.Py() << "," << Tau.Pz() << "," << Tau.E() << ")" << std::endl;
      std::cout << "Sec vertex initial (" << Secvtx_orig.position().x() << "," << Secvtx_orig.position().y() << "," << Secvtx_orig.position().z() 
      	<< ") final (" << Secvtx.position().x() << "," << Secvtx.position().y() << "," << Secvtx.position().z() << ")" << std::endl;

      ///////////////////////////////


      nEvt->Fill(0.5,weight);
      TauMatched->Fill(0.0,weight);
      TauMass->Fill(Tau.M(),weight);
      dTauMass->Fill(Tau.M()-Tau_orig.M(),weight);
      TauPhiChange->Fill(Tau.DeltaPhi(Tau_orig),weight);
      TauThetaChange->Fill(fabs(Tau.Theta()-Tau_orig.Theta()),weight);
      TauEChange->Fill(Tau.E()-Tau_orig.E(),weight);

      for(unsigned int i=0; i<Pions.size();i++){
	PionMass->Fill(Pions.at(i).M(),weight);
	dPionMass->Fill(Pions.at(i).M()-Pions_orig.at(i).M(),weight);
	PionPhiChange->Fill(Pions.at(i).DeltaPhi(Pions_orig.at(i)),weight);
        PionThetaChange->Fill(fabs(Pions.at(i).Theta()-Pions_orig.at(i).Theta()),weight);
        PionEChange->Fill(Pions.at(i).E()-Pions_orig.at(i).E(),weight);
      }

      NuMass->Fill(Nu.M(),weight);
      dNuMass->Fill(Nu.M()-Nu_orig.M(),weight);
      NuPhiChange->Fill(Nu.DeltaPhi(Nu_orig),weight);
      NuThetaChange->Fill(fabs(Nu.Theta()-Nu_orig.Theta()),weight);
      NuEChange->Fill(Nu.Theta()-Nu_orig.Theta(),weight);

      VtxXChange->Fill(Pvtx.position().x()-Pvtx_orig.position().x(),weight);
      VtxYChange->Fill(Pvtx.position().y()-Pvtx_orig.position().y(),weight);
      VtxZChange->Fill(Pvtx.position().z()-Pvtx_orig.position().z(),weight);

      SecVtxXChange->Fill(Secvtx.position().x()-Secvtx_orig.position().x(),weight);
      SecVtxYChange->Fill(Secvtx.position().y()-Secvtx_orig.position().y(),weight);
      SecVtxZChange->Fill(Secvtx.position().z()-Secvtx_orig.position().z(),weight);


      // If Truth is valid run truth comparison
      if(genParticles.isValid()){
	for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	  const reco::GenParticle mytau=(*itr);
	  if(isTruthTauInAcceptance(mytau)){
	    TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	    if(Tau.DeltaR(mc)<TauMatchingDR_){
	      TauMatched->Fill(1.0,weight);
	      TauDecay_CMSSWReco TD;
	      unsigned int jak_id, TauBitMask;
	      TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	      JAKID->Fill(jak_id,weight);
	      std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
	      if(JAKIDtoIndex.count(jak_id)==1){
		unsigned int idx=JAKIDtoIndex.find(jak_id)->second;
		//Tau
		TauMatch_dphi.at(idx)->Fill(Tau.DeltaPhi(mc),weight);
		TauMatch_dtheta.at(idx)->Fill(mc.Theta()-Tau.Theta(),weight);
		TauMatch_e.at(idx)->Fill(mc.E()-Tau.E(),weight);
		//charged hadrons (pi/K)
		for(unsigned int i=0; i<Pions.size();i++){
		  double pidrmin=999;
		  TLorentzVector mcpion;
		  for(unsigned int j=0;j<DecayProd.size();j++){
		    if((fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::pi_plus) ||fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::K_plus)) && DecayProd.at(j)->status()==1){
		      TLorentzVector mcPions_t(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
		      if(mcPions_t.DeltaR(Pions.at(i))<pidrmin){
			mcpion=mcPions_t;
			pidrmin=mcPions_t.DeltaR(Pions.at(i));
		      }
		    }
		  }
		  if(pidrmin<0.1){
		    PionMatch_dphi.at(idx)->Fill(mcpion.DeltaPhi(Pions.at(i)),weight);
		    PionMatch_dtheta.at(idx)->Fill(mcpion.Theta()-Pions.at(i).Theta(),weight);
		    PionMatch_e.at(idx)->Fill(mcpion.E()-Pions.at(i).E(),weight);
		  }
		}
		//nu
		TLorentzVector mcnu;
		bool hasnu=false;
		for(unsigned int i=0;i<DecayProd.size();i++){
		  if(fabs(PdtPdgMini::nu_tau)==fabs(DecayProd.at(i)->pdgId()) && DecayProd.at(i)->status()==1){ 
		    mcnu= TLorentzVector(DecayProd.at(i)->p4().Px(),DecayProd.at(i)->p4().Py(),DecayProd.at(i)->p4().Pz(),DecayProd.at(i)->p4().E()); hasnu=true;
		  }
		}
		if(hasnu){
		  NuMatch_dphi.at(idx)->Fill(mcnu.DeltaPhi(Nu.Phi()),weight);
		  NuMatch_dtheta.at(idx)->Fill(mcnu.Theta()-Nu.Theta(),weight);
		  NuMatch_e.at(idx)->Fill(mcnu.Phi()-Nu.Phi(),weight);
		}
	      }
	    }
	  }
	}
      }
    }
  }


  // Fill truth info for all taus within acceptance
  if(genParticles.isValid()){
    bool hastau=false;
    for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
      const reco::GenParticle mytau=(*itr);
      if(isTruthTauInAcceptance(mytau)){
	hastau=true;
	TauDecay_CMSSWReco TD;
	unsigned int jak_id, TauBitMask;
	TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	JAKIDall->Fill(jak_id,weight);
      }
    }
    if(hastau){
      // nasty hack for Eff - recompute eff every event with a tau
      JAKIDeff->Reset();
      JAKIDeff->getTH1F()->Divide(JAKID->getTH1F(),JAKIDall->getTH1F(), 1., 1., "b");
    }
  }
}



void KinematicTauAnalyzer::beginJob(){
  if(dbe){
    ///Setting the DQM top directories
    dbe->setCurrentFolder("Tau/KinematicFitTau");
    
    // Number of analyzed events
    nEvt = dbe->book1D("nEvt", "n analyzed Events", 1, 0., 1.);                                nEvt->setAxisTitle("Number of Events");
    TauMass           = dbe->book1D("TauMass","M_{Tau}",100,1.7,1.8);                          TauMass->setAxisTitle("M_{Tau} (GeV)");
    PionMass          = dbe->book1D("PionMass","M_{#pi}",100,0.13,0.14);                       PionMass->setAxisTitle("M_{#pi} (GeV)");
    NuMass      = dbe->book1D("NuMass","M_{nu}",100,-0.05,0.05);                               NuMass->setAxisTitle("M_{#nu} (GeV)");

    VtxXChange= dbe->book1D("VtxXChange","Vtx_{X}",100,-0.5,0.5);                              VtxXChange->setAxisTitle("Vtx_{X} (mm)"); 
    VtxYChange= dbe->book1D("VtxYChange","Vtx_{Y}",100,-0.5,0.5);                              VtxYChange->setAxisTitle("Vtx_{Y} (mm)");
    VtxZChange= dbe->book1D("VtxZChange","Vtx_{Z}",100,-0.5,0.5);                              VtxZChange->setAxisTitle("Vtx_{Z} (mm)");
    SecVtxXChange= dbe->book1D("SecVtxXChange","Vtx_{X}^{Sec}",100,-0.02,0.02);                SecVtxXChange->setAxisTitle("Vtx_{X}^{Sec} (mm)");
    SecVtxYChange= dbe->book1D("SecVtxYChange","Vtx_{Y}^{Sec}",100,-0.02,0.02);                SecVtxYChange->setAxisTitle("Vtx_{Y}^{Sec} (mm)");
    SecVtxZChange= dbe->book1D("SecVtxZChange","Vtx_{Z}^{Sec}",100,-0.02,0.02);                SecVtxZChange->setAxisTitle("Vtx_{Z}^{Sec} (mm)");

    TauPhiChange= dbe->book1D("TauPhiChange","#delta#phi_{#tau}",100,-0.5,0.05);               TauPhiChange->setAxisTitle("#delta#phi_{#tau} (rad)");
    TauThetaChange= dbe->book1D("TauThetaChange","#delta#theta_{#tau}",100,-0.05,0.05);        TauThetaChange->setAxisTitle("#delta#theta_{#tau} (rad)");
    TauEChange= dbe->book1D("TauEChange","#deltaE_{#tau}",100,-10,10);                         TauEChange->setAxisTitle("#deltaE_{#tau} (GeV)");
    PionPhiChange= dbe->book1D("PionPhiChange","#delta#phi_{#pi}",100,-0.05,0.05);             PionPhiChange->setAxisTitle("#delta#phi_{#pi} (rad)");
    PionThetaChange= dbe->book1D("PionThetaChange","#delta#theta_{#pi}",100,-0.05,0.05);       PionThetaChange->setAxisTitle("#delta#theta_{#pi} (rad)");
    PionEChange= dbe->book1D("PionEChange","#deltaE_{#pi}",100,-10,10);                        PionEChange->setAxisTitle("#deltaE_{#pi} (GeV)");
    NuPhiChange= dbe->book1D("NuPhiChange","#delta#phi_{#nu}",100,-TMath::Pi(),TMath::Pi());        NuPhiChange->setAxisTitle("#delta#phi_{#nu} (rad)");
    NuThetaChange= dbe->book1D("NuThetaChange","#delta#theta_{#nu}",100,-TMath::Pi(),TMath::Pi());  NuThetaChange->setAxisTitle("#delta#theta_{#nu} (rad)");
    NuEChange= dbe->book1D("NuEChange","#deltaE_{#nu}",100,-10,10);                            NuEChange->setAxisTitle("#deltaE_{#nu} (GeV)");

    dTauMass           = dbe->book1D("dTauMass","M_{Tau}",100,-0.01,0.01);                     dTauMass->setAxisTitle("#deltaM_{Tau} (GeV)");
    dPionMass          = dbe->book1D("dPionMass","M_{#pi}",100,-0.01,0.01);                    dPionMass->setAxisTitle("#deltaM_{#pi} (GeV)");
    dNuMass            = dbe->book1D("dNuMass","M_{nu}",100,-0.05,0.05);                         dNuMass->setAxisTitle("#deltaM_{#nu} (GeV)");

    JAKID =dbe->book1D("JAKID","JAK ID",TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5);            JAKID->setAxisTitle("JAK ID");
    JAKIDall = dbe->book1D("JAKIDall","JAK ID All",TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5); JAKIDall->setAxisTitle("JAK ID");
    JAKIDeff = dbe->book1D("JAKIDeff","JAK ID Eff",TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5); JAKIDeff->setAxisTitle("JAK ID");
    TauMatched = dbe->book1D("TauMatched","TauMatched",2,-0.5,1.5);                                      TauMatched->setAxisTitle("Tau Matched (0=All 1=Matched)");
    for(unsigned int i=0; i<TauDecay::NJAKID;i++){
      if(doJAKID(i)){
	unsigned int idx=TauMatch_dphi.size();
	JAKIDtoIndex.insert(std::pair<unsigned int,unsigned int>(i,idx));
	TString tmp="JAKID";
	tmp+=i;
	TString axis;
	TauMatch_dphi.push_back(dbe->book1D("TauMatchdphi"+tmp,               "d#phi_{"+tmp+"}^{#tau Match}",   100 ,-0.1,0.1));    axis="d#phi_{"+tmp+"}^{#tau Match} (rad)"; TauMatch_dphi.at(idx)->setAxisTitle(axis.Data());
	TauMatch_dtheta.push_back(dbe->book1D("TauMatchdtheta"+tmp,           "d#theta_{"+tmp+"}^{#tau Match}", 100 ,-0.1,0.1));    axis="d#theta_{"+tmp+"}^{#tau Match} (rad)"; TauMatch_dtheta.at(idx)->setAxisTitle(axis.Data());
	TauMatch_e.push_back(dbe->book1D("TauMatche"+tmp,                     "dE_{"+tmp+"}^{#tau Match}",      100 ,-50.0,50.0));  axis="dE_{"+tmp+"}^{#tau Match} (GeV)"; TauMatch_e.at(idx)->setAxisTitle(axis.Data());
	PionMatch_dphi.push_back(dbe->book1D("PionMatchdphi"+tmp,             "d#phi_{"+tmp+"}^{#pi Match}",    100 ,-0.01,0.01));  axis="d#phi_{"+tmp+"}^{#pi Match} (rad)"; PionMatch_dphi.at(idx)->setAxisTitle(axis.Data());
	PionMatch_dtheta.push_back(dbe->book1D("PionMatchdtheta"+tmp,         "d#theta_{"+tmp+"}^{#pi Match}",  100 ,-0.01,0.01));  axis="d#theta_{"+tmp+"}^{#pi Match} (rad)"; PionMatch_dtheta.at(idx)->setAxisTitle(axis.Data());
	PionMatch_e.push_back(dbe->book1D("PionMatche"+tmp,                   "dE_{"+tmp+"}^{#pi Match}",       100 ,-2.0,2.0));    axis="dE_{"+tmp+"}^{#pi Match} (GeV)"; PionMatch_e.at(idx)->setAxisTitle(axis.Data());
	NuMatch_dphi.push_back(dbe->book1D("NuMatchdphi"+tmp,     "d#phi_{"+tmp+"}^{#nu Match}",    100 ,-0.10,0.10));              axis="d#phi_{"+tmp+"}^{#nu Match} (rad)"; NuMatch_dphi.at(idx)->setAxisTitle(axis.Data());
	NuMatch_dtheta.push_back(dbe->book1D("NuMatchdtheta"+tmp, "d#theta_{"+tmp+"}^{#nu Match}",  100 ,-0.20,0.20));              axis="d#theta_{"+tmp+"}^{#nu Match} (rad)"; NuMatch_dtheta.at(idx)->setAxisTitle(axis.Data());
	NuMatch_e.push_back(dbe->book1D("NuMatche"+tmp,           "dE_{"+tmp+"}^{#nu Match}",       100 ,-2.0,2.0));                axis="dE_{"+tmp+"}^{#nu Match} (GeV)"; NuMatch_e.at(idx)->setAxisTitle(axis.Data());
      }
    }
  }
}

void KinematicTauAnalyzer::endJob(){
  //  JAKIDeff->getTH1F()->Divide(JAKID->getTH1F(),JAKIDall->getTH1F(), 1., 1., "b");
}

bool KinematicTauAnalyzer::doJAKID(unsigned int i){
  if(i==TauDecay::JAK_A1_3PI)  return true;
  if(i==TauDecay::JAK_3PIPI0)  return true;
  //if(i==TauDecay::JAK_3PI2PI0) return true;
  //if(i==TauDecay::JAK_3PI3PI0) return true;
  //if(i==TauDecay::JAK_KPIK)    return true;
  //if(i==TauDecay::JAK_KPIPI)   return true;
  return false;
}

bool KinematicTauAnalyzer::isTruthTauInAcceptance(const reco::GenParticle &cand){
  if(fabs(cand.pdgId())!=fabs(PdtPdgMini::tau_minus)) return false;
  if(cand.status()!=2) return false; // require tau that is: 2) a particle after parton showering and ISR/FSR 
  TLorentzVector tau(cand.p4().Px(),cand.p4().Py(),cand.p4().Pz(),cand.p4().E());
  if(tau.Pt()>TauPtMin_ && fabs(tau.Eta())<TauEtaMax_)return true; // require tau within Pt and |eta| acceptance
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauAnalyzer);

