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
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "TVectorT.h"
#include <iostream>

KinematicTauAnalyzer::KinematicTauAnalyzer(const edm::ParameterSet& iConfig):
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  KinematicFitTauTag_(iConfig.getParameter<edm::InputTag>("KinematicFitTauTag")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  TauMatchingDR_( iConfig.getParameter<double>("TauMatchingDR")),
  TauPtMin_( iConfig.getParameter<double>("TauPtMin")),
  TauEtaMax_( iConfig.getParameter<double>("TauEtaMax")),
  JAKID_( iConfig.getParameter< std::vector<int> >("jakid") )
{
  dbe = 0;
  dbe = edm::Service<DQMStore>().operator->();
}

KinematicTauAnalyzer::~KinematicTauAnalyzer(){
}

void KinematicTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  double weight=1.0;
  cnt_++;
  edm::Handle<SelectedKinematicDecayCollection> KinematicFitTaus;
  iEvent.getByLabel(KinematicFitTauTag_,KinematicFitTaus);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(gensrc_, genParticles);

  bool found=false;
  for(SelectedKinematicDecayCollection::const_iterator kinFitTau=KinematicFitTaus->begin();kinFitTau!=KinematicFitTaus->end();kinFitTau++){
    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
      unsigned int npassed=0;
      for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
	if(kinFitTau->discriminators(ambiguity).count(*discr)>0){
	  if(kinFitTau->discriminators(ambiguity).find(*discr)->second) npassed++;
	}
      }
      if(npassed==discriminators_.size()){
	if(found) cntFound_++;
	found=true;
	SelectedKinematicDecay KFTau=(*kinFitTau);
	const TLorentzVector Tau=KFTau.Tau(ambiguity);
	const TLorentzVector a1=KFTau.a1_p4(ambiguity);
	const TLorentzVector a1_initial=KFTau.Initial_a1_p4();
	const TLorentzVector Nu=KFTau.Neutrino(ambiguity);
	const std::vector<TLorentzVector> Pions=KFTau.Pions(ambiguity);

	const TLorentzVector Tau_initial=KFTau.InitialTauGuess(ambiguity);
	const TLorentzVector Nu_initial=KFTau.InitialNeutrinoGuess(ambiguity);
	const std::vector<TLorentzVector> Pions_initial=KFTau.InitialPions();
	
	reco::Vertex Pvtx=KFTau.PrimaryVertexReFitAndRotated();
	reco::Vertex Pvtx_initial=KFTau.InitialPrimaryVertexReFitAndRotated();
	
	reco::Vertex Secvtx=KFTau.SecondaryVertex(ambiguity);
	reco::Vertex Secvtx_initial=KFTau.InitialSecondaryVertex();

	TVector3 FlightDir_initial=KFTau.InitialTauFlghtDirGuess(ambiguity);
	TVector3 FlightDir(Secvtx.position().x()-Pvtx.position().x(),
			   Secvtx.position().y()-Pvtx.position().y(),
			   Secvtx.position().z()-Pvtx.position().z());
	
	nEvt.at(ambiguity)->Fill(0.5,weight);
	Truth_TauMatched.at(ambiguity)->Fill(0.0,weight);
	TauMass.at(ambiguity)->Fill(Tau.M(),weight);
	dTauMass.at(ambiguity)->Fill(Tau.M()-Tau_initial.M(),weight);
	TauPhiChange.at(ambiguity)->Fill(Tau.DeltaPhi(Tau_initial),weight);
	TauThetaChange.at(ambiguity)->Fill(fabs(Tau.Theta()-Tau_initial.Theta()),weight);
	TauEChange.at(ambiguity)->Fill(Tau.E()-Tau_initial.E(),weight);
	
	TauFlightDir.at(ambiguity)->Fill(FlightDir.Angle(Tau.Vect()),weight);
        TauFlightDirInitial.at(ambiguity)->Fill(FlightDir_initial.Angle(Tau_initial.Vect()),weight);

	std::cout << "DQM " << Tau_initial.Px() << " " << Tau_initial.Py() << " " << Tau_initial.Pz() << " " << Tau_initial.E() 
		  << " Vector " << FlightDir_initial.X() << " " << FlightDir_initial.Y() << " " <<FlightDir_initial.Z() << std::endl;

	GFAngle.at(ambiguity)->Fill(a1.Angle(Tau.Vect()),weight);
	GFAngleInitial.at(ambiguity)->Fill(a1_initial.Angle(Tau_initial.Vect()),weight);

	for(unsigned int i=0; i<Pions.size();i++){
	  PionMass.at(ambiguity)->Fill(Pions.at(i).M(),weight);
	  dPionMass.at(ambiguity)->Fill(Pions.at(i).M()-Pions_initial.at(i).M(),weight);
	  PionPhiChange.at(ambiguity)->Fill(Pions.at(i).DeltaPhi(Pions_initial.at(i)),weight);
	  PionThetaChange.at(ambiguity)->Fill(fabs(Pions.at(i).Theta()-Pions_initial.at(i).Theta()),weight);
	  PionEChange.at(ambiguity)->Fill(Pions.at(i).E()-Pions_initial.at(i).E(),weight);
	}
	
	NuMass.at(ambiguity)->Fill(Nu.M(),weight);
	dNuMass.at(ambiguity)->Fill(Nu.M()-Nu_initial.M(),weight);
	NuPhiChange.at(ambiguity)->Fill(Nu.DeltaPhi(Nu_initial),weight);
	NuThetaChange.at(ambiguity)->Fill(fabs(Nu.Theta()-Nu_initial.Theta()),weight);
	NuEChange.at(ambiguity)->Fill(Nu.Theta()-Nu_initial.Theta(),weight);
	
	VtxXChange.at(ambiguity)->Fill(Pvtx.position().x()-Pvtx_initial.position().x(),weight);
	VtxYChange.at(ambiguity)->Fill(Pvtx.position().y()-Pvtx_initial.position().y(),weight);
	VtxZChange.at(ambiguity)->Fill(Pvtx.position().z()-Pvtx_initial.position().z(),weight);

	SecVtxXChange.at(ambiguity)->Fill(Secvtx.position().x()-Secvtx_initial.position().x(),weight);
	SecVtxYChange.at(ambiguity)->Fill(Secvtx.position().y()-Secvtx_initial.position().y(),weight);
	SecVtxZChange.at(ambiguity)->Fill(Secvtx.position().z()-Secvtx_initial.position().z(),weight);
	
	vtxSignPVRotSV.at(ambiguity)->Fill(KFTau.vtxSignPVRotSV(ambiguity),weight);
	vtxSignPVRotPVRed.at(ambiguity)->Fill(KFTau.vtxSignPVRotPVRed(ambiguity),weight);
	a1Mass.at(ambiguity)->Fill(KFTau.a1Mass(ambiguity),weight);
	energyTFraction.at(ambiguity)->Fill(KFTau.energyTFraction(ambiguity),weight);
	iterations.at(ambiguity)->Fill(KFTau.iterations(ambiguity),weight);
        maxiterations.at(ambiguity)->Fill(KFTau.maxiterations(ambiguity),weight);
        chi2.at(ambiguity)->Fill(KFTau.chi2(ambiguity),weight);
        constraints.at(ambiguity)->Fill(KFTau.constraints(ambiguity),weight);
        ndf.at(ambiguity)->Fill(KFTau.ndf(ambiguity),weight);
        csum.at(ambiguity)->Fill(KFTau.csum(ambiguity),weight);
	mincsum.at(ambiguity)->Fill(KFTau.mincsum(ambiguity),weight);
	chi2prob.at(ambiguity)->Fill(KFTau.chi2prob(ambiguity),weight);

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
		Truth_TauMatched.at(ambiguity)->Fill(1.0,weight);
		TauDecay_CMSSWReco TD;
		unsigned int jak_id, TauBitMask;
		TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
		JAKID.at(ambiguity)->Fill(jak_id,weight);
		std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();

		TVector3 TruthPvtx(itr->vx(),itr->vy(),itr->vz());
		TVector3 TruthSvtx;
		for( unsigned int dp=0;dp<DecayProd.size();dp++){
		  if(fabs(DecayProd.at(dp)->pdgId())!=fabs(PdtPdgMini::tau_minus))
		    TruthSvtx=TVector3(DecayProd.at(dp)->vx(),DecayProd.at(dp)->vy(),DecayProd.at(dp)->vz());
		}

		if(JAKIDtoIndex.count(jak_id)==1){
		  unsigned int idx=JAKIDtoIndex.find(jak_id)->second;
                  TMatrixDSym Pvtx_cov(3);
                  Pvtx_cov.ResizeTo(TMatrixDSym(3));
                  for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Pvtx_cov(i,j) = Pvtx.covariance(i,j);//diagonals are squares of sigmas
                  TMatrixDSym Svtx_cov(3);
                  Svtx_cov.ResizeTo(TMatrixDSym(3));
                  for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Svtx_cov(i,j) = Secvtx.covariance(i,j);//diagonals are squares of sigmas
                  TVector3 Pvtx_point(Pvtx.position().x(),Pvtx.position().y(),Pvtx.position().z());
                  TVector3 Svtx_point(Secvtx.position().x(),Secvtx.position().y(),Secvtx.position().z());
                  TMatrixDSym Truth_cov(3);
                  Truth_cov.ResizeTo(TMatrixDSym(3));
                  for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Svtx_cov(i,j) = 0.0;

		  //Tau
		  Truth_TauMatch_dPhi.at(ambiguity).at(idx)->Fill(Tau.DeltaPhi(mc),weight);
		  Truth_TauMatch_dTheta.at(ambiguity).at(idx)->Fill(mc.Theta()-Tau.Theta(),weight);
		  Truth_TauMatch_dE.at(ambiguity).at(idx)->Fill(mc.E()-Tau.E(),weight);
		  Truth_TauMatch_dPt.at(ambiguity).at(idx)->Fill(mc.Pt()-Tau.Pt(),weight);
		  Truth_TauMatch_dPz.at(ambiguity).at(idx)->Fill(mc.Pz()-Tau.Pz(),weight);

                  Truth_TauMatch_dPhiInitial.at(ambiguity).at(idx)->Fill(Tau_initial.DeltaPhi(mc),weight);
                  Truth_TauMatch_dThetaInitial.at(ambiguity).at(idx)->Fill(Tau_initial.Theta()-mc.Theta(),weight);
                  Truth_TauMatch_dEInitial.at(ambiguity).at(idx)->Fill(Tau_initial.E()-mc.E(),weight);

		  TruthVtxXChange.at(ambiguity).at(idx)->Fill(Pvtx.position().x()-TruthPvtx.X(),weight);
		  TruthVtxYChange.at(ambiguity).at(idx)->Fill(Pvtx.position().y()-TruthPvtx.Y(),weight);
		  TruthVtxZChange.at(ambiguity).at(idx)->Fill(Pvtx.position().z()-TruthPvtx.Z(),weight);

		  TruthSecVtxXChange.at(ambiguity).at(idx)->Fill(Secvtx.position().x()-TruthSvtx.X(),weight);
		  TruthSecVtxYChange.at(ambiguity).at(idx)->Fill(Secvtx.position().y()-TruthSvtx.Y(),weight);
		  TruthSecVtxZChange.at(ambiguity).at(idx)->Fill(Secvtx.position().z()-TruthSvtx.Z(),weight);

		  VertexRotation vtxC;
		  TruthPVtxSig.at(ambiguity).at(idx)->Fill(vtxC.vtxDistanceSignificance(Pvtx_point,Pvtx_cov,TruthPvtx,Truth_cov),weight);
		  TruthSecVtxSig.at(ambiguity).at(idx)->Fill(vtxC.vtxDistanceSignificance(Svtx_point,Svtx_cov,TruthSvtx,Truth_cov),weight);

		  TVector3 TruthFlightDir=TruthSvtx-TruthPvtx;
		  TruthTauFlightDir.at(ambiguity).at(idx)->Fill(fabs(Tau.Angle(mc.Vect())),weight);
		  TruthTauFlightDirCheck.at(ambiguity).at(idx)->Fill(fabs(TruthFlightDir.Angle(mc.Vect())),weight);
		  TruthTauFlightDirInitial.at(ambiguity).at(idx)->Fill(fabs(Tau_initial.Angle(mc.Vect())),weight);

                  Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),mc.Pt()-Tau.Pt(),weight);
                  Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),mc.E()-Tau.E(),weight);

		  TLorentzVector mc_a1(0,0,0,0);
		  for(unsigned int j=0;j<DecayProd.size();j++){
		    if(fabs(DecayProd.at(j)->pdgId())>100){
		      TLorentzVector mc_pi(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
		      mc_a1+=mc_pi;
		    }
		  }
		  Truth_TauMatch_dGFAngle.at(ambiguity).at(idx)->Fill(a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFInitialAngle.at(ambiguity).at(idx)->Fill(a1_initial.Angle(Tau_initial.Vect())-mc.Angle(mc_a1.Vect()),weight);
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
		      Truth_PionMatch_dPhi.at(ambiguity).at(idx)->Fill(mcpion.DeltaPhi(Pions.at(i)),weight);
		      Truth_PionMatch_dTheta.at(ambiguity).at(idx)->Fill(mcpion.Theta()-Pions.at(i).Theta(),weight);
		      Truth_PionMatch_dE.at(ambiguity).at(idx)->Fill(mcpion.E()-Pions.at(i).E(),weight);
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
		    Truth_NuMatch_dPhi.at(ambiguity).at(idx)->Fill(mcnu.DeltaPhi(Nu.Phi()),weight);
		    Truth_NuMatch_dTheta.at(ambiguity).at(idx)->Fill(mcnu.Theta()-Nu.Theta(),weight);
		    Truth_NuMatch_dE.at(ambiguity).at(idx)->Fill(mcnu.E()-Nu.E(),weight);
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
  cnt_ = 0;
  cntFound_ = 0;

  if(dbe){
    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
      TString amb="";
      if(ambiguity==SelectedKinematicDecay::ZeroAmbiguitySolution) amb="ZeroAmbiguity";
      if(ambiguity==SelectedKinematicDecay::PlusSolution)          amb="PlusSolution";
      if(ambiguity==SelectedKinematicDecay::MinusSolution)         amb="MinusSolution";
      // Number of analyzed events
      dbe->setCurrentFolder("KinematicFitTau/FitResult");
      nEvt.push_back(dbe->book1D("nEvt"+amb,"n analyzed Events "+amb, 1, 0., 1.));  nEvt.at(ambiguity)->setAxisTitle("Number of Events");
      TauMass.push_back(dbe->book1D("TauMass"+amb,"M_{Tau} "+amb,100,1.7,1.8));      TauMass.at(ambiguity)->setAxisTitle("M_{Tau} (GeV)");
      PionMass.push_back(dbe->book1D("PionMass"+amb,"M_{#pi} "+amb,100,0.13,0.14));  PionMass.at(ambiguity)->setAxisTitle("M_{#pi} (GeV)");
      NuMass.push_back(dbe->book1D("NuMass"+amb,"M_{nu} "+amb,100,-0.05,0.05));      NuMass.at(ambiguity)->setAxisTitle("M_{#nu} (GeV)");
      TauFlightDir.push_back(dbe->book1D("TauFlightDir"+amb,"|#psi_{#tau Dir.,Vtx}^{KF}-#psi_{#tau Dir.,#tau}^{KF}|"+amb,100,0.0,0.5));  TauFlightDir.at(ambiguity)->setAxisTitle("|#psi_{Vtx}^{#tau Dir.,KF}-#psi_{#tau Dir.,#tau}^{KF}| (rad)");
      TauFlightDirInitial.push_back(dbe->book1D("TauFlightDirInitial"+amb,"|#psi_{Vtx}^{Initial}-#psi_{#tau}^{Initial}|"+amb,100,0.0,0.5));  TauFlightDirInitial.at(ambiguity)->setAxisTitle("|#psi_{#tau Dir.,Vtx}^{Initial}-#psi_{#tau Dir.,#tau}^{Initial}| (rad)");

      ////////////////////////////////////////////////////////////////////
      dbe->setCurrentFolder("KinematicFitTau/InitialtoFit");
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

      /////////////////////////////////////////////////////////////////////
      dbe->setCurrentFolder("KinematicFitTau/FitQuality");

      vtxSignPVRotSV.push_back(dbe->book1D("vtxSignPVRotSV"+amb,"vtxSignPVRotSV "+amb,100,0.0,20));  vtxSignPVRotSV.at(ambiguity)->setAxisTitle("#sigma(V_{Prime,Rot},V_{Secondary}) ");
      vtxSignPVRotPVRed.push_back(dbe->book1D("vtxSignPVRotPVRed"+amb,"vtxSignPVRotPVRed "+amb,100,0.0,20.0));  vtxSignPVRotPVRed.at(ambiguity)->setAxisTitle("#sigma(V_{Prime,Rot},V_{Prime}) ");
      a1Mass.push_back(dbe->book1D("a1Mass"+amb,"M_{a1}"+amb,100,0.0,10.0));  a1Mass.at(ambiguity)->setAxisTitle("M_{a1} (GeV)");
      energyTFraction.push_back(dbe->book1D("energyTFraction"+amb,"energyTFraction"+amb,100,0.0,2.0));  energyTFraction.at(ambiguity)->setAxisTitle("E_{a1}/E_{#tau}");
      iterations.push_back(dbe->book1D("iterations"+amb,"iterations"+amb,51,-0.5,50.5));  iterations.at(ambiguity)->setAxisTitle("Number of Iterations");
      maxiterations.push_back(dbe->book1D("maxiterations"+amb,"maxiterations"+amb,51,0.5,50.5));  maxiterations.at(ambiguity)->setAxisTitle("Max Number of Iterations");
      chi2.push_back(dbe->book1D("chi2"+amb,"chi2"+amb,100,0.0,100.0));  chi2.at(ambiguity)->setAxisTitle("#chi^{2}");
      constraints.push_back(dbe->book1D("constraints"+amb,"constraints"+amb,51,-0.5,50.5));  constraints.at(ambiguity)->setAxisTitle("Number of Constraints");
      ndf.push_back(dbe->book1D("ndf"+amb,"ndf"+amb,51,-0.5,50.5));  ndf.at(ambiguity)->setAxisTitle("N.D.F");
      csum.push_back(dbe->book1D("csum"+amb,"csum"+amb,51,-0.5,50.5));  csum.at(ambiguity)->setAxisTitle("csum");
      mincsum.push_back(dbe->book1D("mincsum"+amb,"mincsum"+amb,51,-0.5,50.5));  mincsum.at(ambiguity)->setAxisTitle("mincsum");
      chi2prob.push_back(dbe->book1D("chi2prob"+amb,"chi2prob"+amb,100,0.0,1.0));  chi2prob.at(ambiguity)->setAxisTitle("#chi^{2} Probabilty");
      GFAngleInitial.push_back(dbe->book1D("GFAngleInitial"+amb,"|#theta_{GF}^{Lab,Initial}|"+amb,100,0.0,0.5));  GFAngleInitial.at(ambiguity)->setAxisTitle("|#theta_{GF}^{Lab,Initial}| (rad)");
      GFAngle.push_back(dbe->book1D("GFAngleKF"+amb,"|#theta_{GF}^{Lab,KF}|"+amb,100,0.0,0.5));  GFAngle.at(ambiguity)->setAxisTitle("|#theta_{GF}^{KF,Initial}| (rad)");

      /////////////////////////////////////////////////////////////////////
      // now for Truth
      dbe->setCurrentFolder("KinematicFitTau/Truth");      
      JAKID.push_back(dbe->book1D("JAKID"+amb,"JAK ID "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));            JAKID.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDall.push_back(dbe->book1D("JAKIDall"+amb,"JAK ID All "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5)); JAKIDall.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDeff.push_back(dbe->book1D("JAKIDeff"+amb,"JAK ID Eff "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5)); JAKIDeff.at(ambiguity)->setAxisTitle("JAK ID");
      Truth_TauMatched.push_back(dbe->book1D("Truth_TauMatched"+amb,"Truth_TauMatched "+amb,2,-0.5,1.5));                                      Truth_TauMatched.at(ambiguity)->setAxisTitle("Tau Matched (0=All 1=Matched)");

      Truth_TauMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dE.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPhiInitial.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dThetaInitial.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEInitial.push_back(std::vector<MonitorElement*>());


      Truth_PionMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dE.push_back(std::vector<MonitorElement*>());
      Truth_NuMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_NuMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_NuMatch_dE.push_back(std::vector<MonitorElement*>());

      TruthVtxXChange.push_back(std::vector<MonitorElement*>());
      TruthVtxYChange.push_back(std::vector<MonitorElement*>());
      TruthVtxZChange.push_back(std::vector<MonitorElement*>());
      TruthSecVtxXChange.push_back(std::vector<MonitorElement*>());
      TruthSecVtxYChange.push_back(std::vector<MonitorElement*>());
      TruthSecVtxZChange.push_back(std::vector<MonitorElement*>());
      TruthPVtxSig.push_back(std::vector<MonitorElement*>());
      TruthSecVtxSig.push_back(std::vector<MonitorElement*>());
      TruthTauFlightDir.push_back(std::vector<MonitorElement*>());
      TruthTauFlightDirCheck.push_back(std::vector<MonitorElement*>());
      TruthTauFlightDirInitial.push_back(std::vector<MonitorElement*>());

      Truth_TauMatch_dPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPz.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPtvsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAngle.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFInitialAngle.push_back(std::vector<MonitorElement*>());

      for(unsigned int i=0; i<TauDecay::NJAKID;i++){
	if(doJAKID(i)){
	  unsigned int idx=Truth_TauMatch_dPhi.at(ambiguity).size();
	  if(JAKIDtoIndex.count(i)==0) JAKIDtoIndex.insert(std::pair<unsigned int,unsigned int>(i,idx));
	  TString tmp="JAKID";
	  tmp+=i;
	  TString axis;
	  Truth_TauMatch_dPhi.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdPhi"+tmp+amb,"d#phi_{"+tmp+"}^{#tau Match} "+amb,100 ,-0.1,0.1));          axis="d#phi_{"+tmp+"}^{#tau Match} (rad)"; Truth_TauMatch_dPhi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_TauMatch_dTheta.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdTheta"+tmp+amb,"d#theta_{"+tmp+"}^{#tau Match} "+amb,100 ,-0.1,0.1));    axis="d#theta_{"+tmp+"}^{#tau Match} (rad)"; Truth_TauMatch_dTheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_TauMatch_dE.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));                 axis="dE_{"+tmp+"}^{#tau Match} (GeV)"; Truth_TauMatch_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_PionMatch_dPhi.at(ambiguity).push_back(dbe->book1D("TruthPionMatchdPhi"+tmp+amb,"d#phi_{"+tmp+"}^{#pi Match} "+amb,100 ,-0.01,0.01));       axis="d#phi_{"+tmp+"}^{#pi Match} (rad)"; Truth_PionMatch_dPhi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_PionMatch_dTheta.at(ambiguity).push_back(dbe->book1D("TruthPionMatchdTheta"+tmp+amb,"d#theta_{"+tmp+"}^{#pi Match} "+amb,100 ,-0.01,0.01)); axis="d#theta_{"+tmp+"}^{#pi Match} (rad)"; Truth_PionMatch_dTheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_PionMatch_dE.at(ambiguity).push_back(dbe->book1D("TruthPionMatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{#pi Match} "+amb,100 ,-2.0,2.0));                  axis="dE_{"+tmp+"}^{#pi Match} (GeV)"; Truth_PionMatch_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_NuMatch_dPhi.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdPhi"+tmp+amb,"d#phi_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.10,0.10));           axis="d#phi_{"+tmp+"}^{#nu Match} (rad)"; Truth_NuMatch_dPhi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_NuMatch_dTheta.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdTheta"+tmp+amb, "d#theta_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.20,0.20));    axis="d#theta_{"+tmp+"}^{#nu Match} (rad)"; Truth_NuMatch_dTheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_NuMatch_dE.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{#nu Match} "+amb,100 ,-2.0,2.0));                      axis="dE_{"+tmp+"}^{#nu Match} (GeV)"; Truth_NuMatch_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());


	  TruthVtxXChange.at(ambiguity).push_back(dbe->book1D("TruthPVtxXChange"+tmp+amb,"Vtx_{x,"+tmp+"}^{Prime,KF}-Vtx_{x}^{Prime,Truth} "+amb,100 ,-0.1,0.1)); axis="Vtx_{x,"+tmp+"}^{Prime,KF}-Vtx_{x}^{Prime,Truth} (mm)"; TruthVtxXChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthVtxYChange.at(ambiguity).push_back(dbe->book1D("TruthPVtxYChange"+tmp+amb,"Vtx_{y,"+tmp+"}^{Prime,KF}-Vtx_{y}^{Prime,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{y,"+tmp+"}^{Prime,KF}-Vtx_{y}^{Prime,Truth} (mm)"; TruthVtxYChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthVtxZChange.at(ambiguity).push_back(dbe->book1D("TruthPVtxZChange"+tmp+amb,"Vtx_{z,"+tmp+"}^{Prime,KF}-Vtx_{z}^{Prime,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{z,"+tmp+"}^{Prime,KF}-Vtx_{z}^{Prime,Truth} (mm)"; TruthVtxZChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthSecVtxXChange.at(ambiguity).push_back(dbe->book1D("TruthSecVtxXChange"+tmp+amb,"Vtx_{x,"+tmp+"}^{Sec.,KF}-Vtx_{x}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{x,"+tmp+"}^{Sec.,KF}-Vtx_{x}^{Sec.,Truth} (mm)"; TruthSecVtxXChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxYChange.at(ambiguity).push_back(dbe->book1D("TruthSecVtxYChange"+tmp+amb,"Vtx_{y,"+tmp+"}^{Sec,KF}-Vtx_{y}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{y,"+tmp+"}^{Sec.,KF}-Vtx_{y}^{Sec.,Truth} (mm)"; TruthSecVtxYChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxZChange.at(ambiguity).push_back(dbe->book1D("TruthSecVtxZChange"+tmp+amb,"Vtx_{z,"+tmp+"}^{Sec,KF}-Vtx_{z}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{z,"+tmp+"}^{Sec.,KF}-Vtx_{z}^{Sec.,Truth} (mm)"; TruthSecVtxZChange.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

	  TruthPVtxSig.at(ambiguity).push_back(dbe->book1D("TruthPVtxSig"+tmp+amb,"#sigma_{Prime Vtx,"+tmp+"}^{Truth} "+amb,100 ,0.0,10.0)); axis="#sigma_{Prime Vtx,"+tmp+"}^{Truth}"; TruthPVtxSig.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthSecVtxSig.at(ambiguity).push_back(dbe->book1D("TruthSecVtxSig"+tmp+amb,"#sigma_{Sec. Vtx"+tmp+"}^{Truth} "+amb,100 ,0.0,10.0));axis="#sigma_{Sec. Vtx"+tmp+"}^{Truth} "; TruthSecVtxSig.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthTauFlightDir.at(ambiguity).push_back(dbe->book1D("TruthTauFlightDir"+tmp+amb,"|d#psi_{#tau Dir.,#tau}^{KF}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| "+amb,100 ,0.0,0.1));axis="|d#psi_{#tau Dir.,#tau}^{KF}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| (rad)"; TruthTauFlightDir.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthTauFlightDirCheck.at(ambiguity).push_back(dbe->book1D("TruthTauFlightDirCheck"+tmp+amb,"|#psi_{#tau Dir.,Vtx}^{Truth}-#psi_{#tau Dir.,#tau}^{Truth}| "+amb,100 ,0.0,0.01));axis="|#psi_{#tau Dir.,Vtx}^{Truth}-#psi_{#tau Dir.,#tau}^{Truth}| (rad)"; TruthTauFlightDirCheck.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthTauFlightDirInitial.at(ambiguity).push_back(dbe->book1D("TruthTauFlightDirInitial"+tmp+amb,"|d#psi_{#tau Dir.,#tau}^{Initial}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| "+amb,100 ,0.0,0.1));axis="|d#psi_{#tau Dir.,#tau}^{Initial}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| (rad)"; TruthTauFlightDirInitial.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

	  Truth_TauMatch_dPt.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdPt"+tmp+amb,"dP_{t,"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));   axis="dP_{t,"+tmp+"}^{#tau Match} (GeV)"; Truth_TauMatch_dPt.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_TauMatch_dPz.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdPz"+tmp+amb,"dP_{z,"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));   axis="dP_{z,"+tmp+"}^{#tau Match} (GeV)"; Truth_TauMatch_dPz.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

	  Truth_TauMatch_dPtvsL.at(ambiguity).push_back(dbe->book2D("TruthTauMatchdPtvsL"+tmp+amb,"Truth Tau-Matchd P_{t} vs L( "+tmp+amb+")",10,0.0,5.0,100,-50,50)); Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->setAxisTitle("L (mm)",1);
	  Truth_TauMatch_dEvsL.at(ambiguity).push_back(dbe->book2D("TruthTauMatchdEvsL"+tmp+amb,"Truth Tau Matchd E vs L "+tmp+amb+")",100,0.0,5.0,100,-50,50)); Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->setAxisTitle("L (mm)",1);

	  Truth_TauMatch_dGFAnglevsL.at(ambiguity).push_back(dbe->book2D("TruthTauMatchddGFAnglevsL"+tmp+amb,"Truth Tau-Matchd d#theta_{GF}^{Lab} vs L( "+tmp+amb+")",10,0.0,5.0,64,-0.1,0.1)); Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (GeV)",2); Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->setAxisTitle("L (mm)",1);

	  Truth_TauMatch_dGFAngle.at(ambiguity).push_back(dbe->book1D("TruthTauMatchddGFAngle"+tmp+amb,"Truth Tau-Matchd d#theta_{GF}^{Lab} "+tmp+amb+")",100,-0.1,0.1)); Truth_TauMatch_dGFAngle.at(ambiguity).at(idx)->setAxisTitle("d#theta_{GF,#tau} (GeV)"); 
	  Truth_TauMatch_dGFInitialAngle.at(ambiguity).push_back(dbe->book1D("TruthTauMatchddGFInitialAngle"+tmp+amb,"Truth Tau-Matchd d#theta_{GF}^{Initial,Lab} "+tmp+amb+")",100,-0.1,0.1)); Truth_TauMatch_dGFAngle.at(ambiguity).at(idx)->setAxisTitle("d#theta_{GF,#tau}^{Initial} (GeV)");


	  Truth_TauMatch_dPhiInitial.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdPhiInitial"+tmp+amb,"d#phi_{"+tmp+"}^{#tau Match,Initial} "+amb,100 ,-0.1,0.1));          axis="d#phi_{"+tmp+"}^{#tau Match,Initial} (rad)"; Truth_TauMatch_dPhiInitial.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          Truth_TauMatch_dThetaInitial.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdThetaInitial"+tmp+amb,"d#theta_{"+tmp+"}^{#tau Match,Initial} "+amb,100 ,-0.1,0.1));    axis="d#theta_{"+tmp+"}^{#tau Match,Initial} (rad)"; Truth_TauMatch_dThetaInitial.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          Truth_TauMatch_dEInitial.at(ambiguity).push_back(dbe->book1D("TruthTauMatchdEnergyInitial"+tmp+amb,"dE_{"+tmp+"}^{#tau Match,Initial} "+amb,100 ,-50.0,50.0));                 axis="dE_{"+tmp+"}^{#tau Match,Initial} (GeV)"; Truth_TauMatch_dEInitial.at(ambiguity).at(idx)->setAxisTitle(axis.Data());


	} 
      }
    }
  }
}

void KinematicTauAnalyzer::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  std::cout << "Tau_JAKID_Filter" <<"--> [Tau_JAKID_Filter]  Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%" << std::endl;
}

bool KinematicTauAnalyzer::doJAKID(unsigned int i){
  for(unsigned int j=0;j<JAKID_.size();j++){
    if(((unsigned int)JAKID_.at(j))==i)return true;
  }
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

