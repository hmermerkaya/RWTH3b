
/**
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */
#include "RecoTauTag/KinematicTau/interface/KinematicTauAnalyzer.h"
// object includes
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
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
#include <iostream>

KinematicTauAnalyzer::KinematicTauAnalyzer(const edm::ParameterSet& iConfig):
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  KinematicFitTauTag_(iConfig.getParameter<edm::InputTag>("KinematicFitTauTag")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  TauMatchingDR_( iConfig.getParameter<double>("TauMatchingDR")),
  TauPtMin_( iConfig.getParameter<double>("TauPtMin")),
  TauEtaMax_( iConfig.getParameter<double>("TauEtaMax"))
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
    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
      unsigned int npassed=0;
      for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
	if(kinFitTau->discriminators(ambiguity).count(*discr)>0){
	  if(kinFitTau->discriminators(ambiguity).find(*discr)->second) npassed++;
	}
      }
      if(npassed==discriminators_.size()){
	SelectedKinematicDecay KFTau=(*kinFitTau);
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
	
	/*	const SelectedKinematicParticleCollection& Particles =KFTau.particles();
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
	*/
	///////////////////////////////
	
	
	nEvt.at(ambiguity)->Fill(0.5,weight);
	TauMatched.at(ambiguity)->Fill(0.0,weight);
	TauMass.at(ambiguity)->Fill(Tau.M(),weight);
	dTauMass.at(ambiguity)->Fill(Tau.M()-Tau_orig.M(),weight);
	TauPhiChange.at(ambiguity)->Fill(Tau.DeltaPhi(Tau_orig),weight);
	TauThetaChange.at(ambiguity)->Fill(fabs(Tau.Theta()-Tau_orig.Theta()),weight);
	TauEChange.at(ambiguity)->Fill(Tau.E()-Tau_orig.E(),weight);
	
	for(unsigned int i=0; i<Pions.size();i++){
	  PionMass.at(ambiguity)->Fill(Pions.at(i).M(),weight);
	  dPionMass.at(ambiguity)->Fill(Pions.at(i).M()-Pions_orig.at(i).M(),weight);
	  PionPhiChange.at(ambiguity)->Fill(Pions.at(i).DeltaPhi(Pions_orig.at(i)),weight);
	  PionThetaChange.at(ambiguity)->Fill(fabs(Pions.at(i).Theta()-Pions_orig.at(i).Theta()),weight);
	  PionEChange.at(ambiguity)->Fill(Pions.at(i).E()-Pions_orig.at(i).E(),weight);
	}
	
	NuMass.at(ambiguity)->Fill(Nu.M(),weight);
	dNuMass.at(ambiguity)->Fill(Nu.M()-Nu_orig.M(),weight);
	NuPhiChange.at(ambiguity)->Fill(Nu.DeltaPhi(Nu_orig),weight);
	NuThetaChange.at(ambiguity)->Fill(fabs(Nu.Theta()-Nu_orig.Theta()),weight);
	NuEChange.at(ambiguity)->Fill(Nu.Theta()-Nu_orig.Theta(),weight);
	
	VtxXChange.at(ambiguity)->Fill(Pvtx.position().x()-Pvtx_orig.position().x(),weight);
	VtxYChange.at(ambiguity)->Fill(Pvtx.position().y()-Pvtx_orig.position().y(),weight);
	VtxZChange.at(ambiguity)->Fill(Pvtx.position().z()-Pvtx_orig.position().z(),weight);

	SecVtxXChange.at(ambiguity)->Fill(Secvtx.position().x()-Secvtx_orig.position().x(),weight);
	SecVtxYChange.at(ambiguity)->Fill(Secvtx.position().y()-Secvtx_orig.position().y(),weight);
	SecVtxZChange.at(ambiguity)->Fill(Secvtx.position().z()-Secvtx_orig.position().z(),weight);

	
	// If Truth is valid run truth comparison
	if(genParticles.isValid()){
	  for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	    const reco::GenParticle mytau=(*itr);
	    if(isTruthTauInAcceptance(mytau)){
	      TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	      if(Tau.DeltaR(mc)<TauMatchingDR_){
		if(ambiguity==SelectedKinematicDecay::PlusSolution || ambiguity==SelectedKinematicDecay::MinusSolution){
		  TLorentzVector Tau_plus=KFTau.Tau(SelectedKinematicDecay::PlusSolution);
		  TLorentzVector Tau_minus=KFTau.Tau(SelectedKinematicDecay::MinusSolution);
		  if(fabs(mc.E()-Tau_plus.E())<fabs(mc.E()-Tau_minus.E()) && Tau_plus.E()!=Tau.E())  continue;
		  if(fabs(mc.E()-Tau_plus.E())>fabs(mc.E()-Tau_minus.E()) && Tau_minus.E()!=Tau.E()) continue;
		}
		TauMatched.at(ambiguity)->Fill(1.0,weight);
		TauDecay_CMSSWReco TD;
		unsigned int jak_id, TauBitMask;
		TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
		JAKID.at(ambiguity)->Fill(jak_id,weight);
		std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
		if(JAKIDtoIndex.count(jak_id)==1){
		  unsigned int idx=JAKIDtoIndex.find(jak_id)->second;
		  //Tau
		  TauMatch_dphi.at(ambiguity).at(idx)->Fill(Tau.DeltaPhi(mc),weight);
		  TauMatch_dtheta.at(ambiguity).at(idx)->Fill(mc.Theta()-Tau.Theta(),weight);
		  TauMatch_e.at(ambiguity).at(idx)->Fill(mc.E()-Tau.E(),weight);
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
		      PionMatch_dphi.at(ambiguity).at(idx)->Fill(mcpion.DeltaPhi(Pions.at(i)),weight);
		      PionMatch_dtheta.at(ambiguity).at(idx)->Fill(mcpion.Theta()-Pions.at(i).Theta(),weight);
		      PionMatch_e.at(ambiguity).at(idx)->Fill(mcpion.E()-Pions.at(i).E(),weight);
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
		    NuMatch_dphi.at(ambiguity).at(idx)->Fill(mcnu.DeltaPhi(Nu.Phi()),weight);
		    NuMatch_dtheta.at(ambiguity).at(idx)->Fill(mcnu.Theta()-Nu.Theta(),weight);
		    NuMatch_e.at(ambiguity).at(idx)->Fill(mcnu.Phi()-Nu.Phi(),weight);
		  }
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
    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
      bool hastau=false;
      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	const reco::GenParticle mytau=(*itr);
	if(isTruthTauInAcceptance(mytau)){
	  hastau=true;
	  TauDecay_CMSSWReco TD;
	  unsigned int jak_id, TauBitMask;
	  TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	  JAKIDall.at(ambiguity)->Fill(jak_id,weight);
	}
      }
      if(hastau){
	// nasty hack for Eff - recompute eff every event with a tau
	JAKIDeff.at(ambiguity)->Reset();
	JAKIDeff.at(ambiguity)->getTH1F()->Divide(JAKID.at(ambiguity)->getTH1F(),JAKIDall.at(ambiguity)->getTH1F(), 1., 1., "b");
      }
    }
  }
}



void KinematicTauAnalyzer::beginJob(){
  if(dbe){
    ///Setting the DQM top directories
    dbe->setCurrentFolder("Tau/KinematicFitTau");
    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
      TString amb="";
      if(ambiguity==SelectedKinematicDecay::ZeroAmbiguitySolution) amb="ZeroAmbiguity";
      if(ambiguity==SelectedKinematicDecay::PlusSolution)          amb="PlusSolution";
      if(ambiguity==SelectedKinematicDecay::MinusSolution)         amb="MinusSolution";
      // Number of analyzed events
      nEvt.push_back(dbe->book1D("nEvt"+amb,"n analyzed Events "+amb, 1, 0., 1.));  nEvt.at(ambiguity)->setAxisTitle("Number of Events");
      TauMass.push_back(dbe->book1D("TauMass"+amb,"M_{Tau} "+amb,100,1.7,1.8));      TauMass.at(ambiguity)->setAxisTitle("M_{Tau} (GeV)");
      PionMass.push_back(dbe->book1D("PionMass"+amb,"M_{#pi} "+amb,100,0.13,0.14));  PionMass.at(ambiguity)->setAxisTitle("M_{#pi} (GeV)");
      NuMass.push_back(dbe->book1D("NuMass"+amb,"M_{nu} "+amb,100,-0.05,0.05));      NuMass.at(ambiguity)->setAxisTitle("M_{#nu} (GeV)");
      VtxXChange.push_back(dbe->book1D("VtxXChange"+amb,"Vtx_{X} "+amb,100,-0.5,0.5));       VtxXChange.at(ambiguity)->setAxisTitle("Vtx_{X} (mm)"); 
      VtxYChange.push_back(dbe->book1D("VtxYChange"+amb,"Vtx_{Y} "+amb,100,-0.5,0.5));       VtxYChange.at(ambiguity)->setAxisTitle("Vtx_{Y} (mm)");
      VtxZChange.push_back(dbe->book1D("VtxZChange"+amb,"Vtx_{Z} "+amb,100,-0.5,0.5));       VtxZChange.at(ambiguity)->setAxisTitle("Vtx_{Z} (mm)");
      SecVtxXChange.push_back(dbe->book1D("SecVtxXChange"+amb,"Vtx_{X}^{Sec} "+amb,100,-0.02,0.02)); SecVtxXChange.at(ambiguity)->setAxisTitle("Vtx_{X}^{Sec} (mm)");
      SecVtxYChange.push_back(dbe->book1D("SecVtxYChange"+amb,"Vtx_{Y}^{Sec} "+amb,100,-0.02,0.02)); SecVtxYChange.at(ambiguity)->setAxisTitle("Vtx_{Y}^{Sec} (mm)");
      SecVtxZChange.push_back(dbe->book1D("SecVtxZChange"+amb,"Vtx_{Z}^{Sec} "+amb,100,-0.02,0.02)); SecVtxZChange.at(ambiguity)->setAxisTitle("Vtx_{Z}^{Sec} (mm)");
      
      TauPhiChange.push_back(dbe->book1D("TauPhiChange"+amb,"#delta#phi_{#tau} "+amb,100,-0.5,0.05));         TauPhiChange.at(ambiguity)->setAxisTitle("#delta#phi_{#tau} (rad)");
      TauThetaChange.push_back(dbe->book1D("TauThetaChange"+amb,"#delta#theta_{#tau} "+amb,100,-0.05,0.05));  TauThetaChange.at(ambiguity)->setAxisTitle("#delta#theta_{#tau} (rad)");
      TauEChange.push_back(dbe->book1D("TauEChange"+amb,"#deltaE_{#tau} "+amb,100,-10,10));                   TauEChange.at(ambiguity)->setAxisTitle("#deltaE_{#tau} (GeV)");
      PionPhiChange.push_back(dbe->book1D("PionPhiChange"+amb,"#delta#phi_{#pi} "+amb,100,-0.05,0.05));       PionPhiChange.at(ambiguity)->setAxisTitle("#delta#phi_{#pi} (rad)");
      PionThetaChange.push_back(dbe->book1D("PionThetaChange"+amb,"#delta#theta_{#pi} "+amb,100,-0.05,0.05)); PionThetaChange.at(ambiguity)->setAxisTitle("#delta#theta_{#pi} (rad)");
      PionEChange.push_back(dbe->book1D("PionEChange"+amb,"#deltaE_{#pi} "+amb,100,-10,10));                  PionEChange.at(ambiguity)->setAxisTitle("#deltaE_{#pi} (GeV)");
      NuPhiChange.push_back(dbe->book1D("NuPhiChange"+amb,"#delta#phi_{#nu} "+amb,100,-TMath::Pi(),TMath::Pi()));        NuPhiChange.at(ambiguity)->setAxisTitle("#delta#phi_{#nu} (rad)");
      NuThetaChange.push_back(dbe->book1D("NuThetaChange"+amb,"#delta#theta_{#nu} "+amb,100,-TMath::Pi(),TMath::Pi()));  NuThetaChange.at(ambiguity)->setAxisTitle("#delta#theta_{#nu} (rad)");
      NuEChange.push_back(dbe->book1D("NuEChange"+amb,"#deltaE_{#nu} "+amb,100,-10,10));                            NuEChange.at(ambiguity)->setAxisTitle("#deltaE_{#nu} (GeV)");

      dTauMass.push_back(dbe->book1D("dTauMass"+amb,"M_{Tau} "+amb,100,-0.01,0.01));                     dTauMass.at(ambiguity)->setAxisTitle("#deltaM_{Tau} (GeV)");
      dPionMass.push_back(dbe->book1D("dPionMass"+amb,"M_{#pi} "+amb,100,-0.01,0.01));                    dPionMass.at(ambiguity)->setAxisTitle("#deltaM_{#pi} (GeV)");
      dNuMass.push_back(dbe->book1D("dNuMass"+amb,"M_{nu} "+amb,100,-0.05,0.05));                         dNuMass.at(ambiguity)->setAxisTitle("#deltaM_{#nu} (GeV)");
      
      JAKID.push_back(dbe->book1D("JAKID"+amb,"JAK ID "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));            JAKID.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDall.push_back(dbe->book1D("JAKIDall"+amb,"JAK ID All "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5)); JAKIDall.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDeff.push_back(dbe->book1D("JAKIDeff"+amb,"JAK ID Eff "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5)); JAKIDeff.at(ambiguity)->setAxisTitle("JAK ID");
      TauMatched.push_back(dbe->book1D("TauMatched"+amb,"TauMatched "+amb,2,-0.5,1.5));                                      TauMatched.at(ambiguity)->setAxisTitle("Tau Matched (0=All 1=Matched)");
      //////////////////////////////////////////////////////////////////
      // now for the truth
      TauMatch_dphi.push_back(std::vector<MonitorElement*>());
      TauMatch_dtheta.push_back(std::vector<MonitorElement*>());
      TauMatch_e.push_back(std::vector<MonitorElement*>());
      PionMatch_dphi.push_back(std::vector<MonitorElement*>());
      PionMatch_dtheta.push_back(std::vector<MonitorElement*>());
      PionMatch_e.push_back(std::vector<MonitorElement*>());
      NuMatch_dphi.push_back(std::vector<MonitorElement*>());
      NuMatch_dtheta.push_back(std::vector<MonitorElement*>());
      NuMatch_e.push_back(std::vector<MonitorElement*>());
      for(unsigned int i=0; i<TauDecay::NJAKID;i++){
	if(doJAKID(i)){
	  unsigned int idx=TauMatch_dphi.at(ambiguity).size();
	  if(JAKIDtoIndex.count(i)==0) JAKIDtoIndex.insert(std::pair<unsigned int,unsigned int>(i,idx));
	  TString tmp="JAKID";
	  tmp+=i;
	  TString axis;
	  TauMatch_dphi.at(ambiguity).push_back(dbe->book1D("TauMatchdphi"+tmp+amb,"d#phi_{"+tmp+"}^{#tau Match} "+amb,100 ,-0.1,0.1));          axis="d#phi_{"+tmp+"}^{#tau Match} (rad)"; TauMatch_dphi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TauMatch_dtheta.at(ambiguity).push_back(dbe->book1D("TauMatchdtheta"+tmp+amb,"d#theta_{"+tmp+"}^{#tau Match} "+amb,100 ,-0.1,0.1));    axis="d#theta_{"+tmp+"}^{#tau Match} (rad)"; TauMatch_dtheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TauMatch_e.at(ambiguity).push_back(dbe->book1D("TauMatche"+tmp+amb,"dE_{"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));                 axis="dE_{"+tmp+"}^{#tau Match} (GeV)"; TauMatch_e.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  PionMatch_dphi.at(ambiguity).push_back(dbe->book1D("PionMatchdphi"+tmp+amb,"d#phi_{"+tmp+"}^{#pi Match} "+amb,100 ,-0.01,0.01));       axis="d#phi_{"+tmp+"}^{#pi Match} (rad)"; PionMatch_dphi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  PionMatch_dtheta.at(ambiguity).push_back(dbe->book1D("PionMatchdtheta"+tmp+amb,"d#theta_{"+tmp+"}^{#pi Match} "+amb,100 ,-0.01,0.01)); axis="d#theta_{"+tmp+"}^{#pi Match} (rad)"; PionMatch_dtheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  PionMatch_e.at(ambiguity).push_back(dbe->book1D("PionMatche"+tmp+amb,"dE_{"+tmp+"}^{#pi Match} "+amb,100 ,-2.0,2.0));                  axis="dE_{"+tmp+"}^{#pi Match} (GeV)"; PionMatch_e.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  NuMatch_dphi.at(ambiguity).push_back(dbe->book1D("NuMatchdphi"+tmp+amb,"d#phi_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.10,0.10));           axis="d#phi_{"+tmp+"}^{#nu Match} (rad)"; NuMatch_dphi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  NuMatch_dtheta.at(ambiguity).push_back(dbe->book1D("NuMatchdtheta"+tmp+amb, "d#theta_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.20,0.20));    axis="d#theta_{"+tmp+"}^{#nu Match} (rad)"; NuMatch_dtheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  NuMatch_e.at(ambiguity).push_back(dbe->book1D("NuMatche"+tmp+amb,"dE_{"+tmp+"}^{#nu Match} "+amb,100 ,-2.0,2.0));                      axis="dE_{"+tmp+"}^{#nu Match} (GeV)"; NuMatch_e.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	}
      }
    }
  }
}

void KinematicTauAnalyzer::endJob(){
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

