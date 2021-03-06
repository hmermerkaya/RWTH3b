/*
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010

 Modifed by Ian M. Nugent
 */
#include "RecoTauTag/KinematicTau/interface/KinematicTauAnalyzer.h"
// object includes
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
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
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "TMath.h"
#include "TMatrixTSym.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "RecoTauTag/KinematicTau/interface/Tau_JAKID_Filter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
 
KinematicTauAnalyzer::KinematicTauAnalyzer(const edm::ParameterSet& iConfig):
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  KinematicFitTauTag_(iConfig.getParameter<edm::InputTag>("KinematicFitTauTag")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  tauType_( iConfig.getUntrackedParameter<std::string>("tauType", "hps") ),
  TauMatchingDR_( iConfig.getParameter<double>("TauMatchingDR")),
  TauPtMin_( iConfig.getParameter<double>("TauPtMin")),
  TauEtaMax_( iConfig.getParameter<double>("TauEtaMax")),
  JAKID_( iConfig.getParameter< std::vector<int> >("jakid") ),
  doFakeRate_( iConfig.getParameter<bool>("doFakeRate")),
  doDQM_( iConfig.getParameter<bool>("doDQM"))
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

  if(doFakeRate_){
    edm::Handle<reco::PFTauCollection> inputCollection;
    iEvent.getByLabel(tauType_+"PFTauProducer", inputCollection);

    for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
      for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
	reco::PFTauRef thePFTau(inputCollection, iPFTau);
	if(fabs(thePFTau->eta())<TauEtaMax_ && thePFTau->pt()>TauPtMin_){
	FakeRate_eta_All.at(ambiguity)->Fill(thePFTau->eta(),weight);
	FakeRate_pt_All.at(ambiguity)->Fill(thePFTau->pt(),weight);
	for(SelectedKinematicDecayCollection::const_iterator kinFitTau=KinematicFitTaus->begin();kinFitTau!=KinematicFitTaus->end();kinFitTau++){
	  SelectedKinematicDecay KFTau=(*kinFitTau);
	  if(pow(KFTau.PFTauRef()->eta()-thePFTau->eta(),2.0)+pow(fmod(fabs(KFTau.PFTauRef()->phi()-thePFTau->phi()),2*TMath::Pi()),2.0)<0.2){
	    FakeRate_eta_isFit.at(ambiguity)->Fill(thePFTau->eta(),weight);
	    FakeRate_pt_isFit.at(ambiguity)->Fill(thePFTau->pt(),weight);
	     bool passed=true;
	    for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
	      if(kinFitTau->discriminators(ambiguity).count(*discr)>0){
		if(!kinFitTau->discriminators(ambiguity).find(*discr)->second) passed=false;
	      }
	    }
	    if(passed){
	      FakeRate_eta.at(ambiguity)->Fill(thePFTau->eta(),weight);
	      FakeRate_pt.at(ambiguity)->Fill(thePFTau->pt(),weight);
	      const std::vector<TLorentzVector> Pions=KFTau.Pions(ambiguity);
	      if(Pions.size()==3){
		std::vector<int> Charges=KFTau.Pions_Charge();
		int charge=0;
		for(unsigned int i=0;i<Pions.size();i++){
		  charge+=Charges.at(i);
		}
		TLorentzVector LV1, LV2, LV3;
		bool has1=false;;
		for(unsigned int i=0;i<Pions.size();i++){
		  if(has1==false && Charges.at(i)==charge){ LV1=Pions.at(i); has1=true;}
		  else if(Charges.at(i)==charge){ LV2=Pions.at(i);}
		  else{ LV3=Pions.at(i); }
		}		
		TLorentzVector LV13=LV1+LV3;
		TLorentzVector LV23=LV2+LV3;
                TLorentzVector LV12=LV1+LV2;
		if(LV13.M()>LV23.M()){
		  DalitzFake.at(ambiguity)->Fill(LV13.M(),LV23.M(),weight);
		  InvMass12Fake.at(ambiguity)->Fill(LV12.M(),weight);
		  InvMass13Fake.at(ambiguity)->Fill(LV13.M(),weight);
		  InvMass23Fake.at(ambiguity)->Fill(LV23.M(),weight);
		}
		else{
		  DalitzFake.at(ambiguity)->Fill(LV23.M(),LV13.M(),weight);
                  InvMass12Fake.at(ambiguity)->Fill(LV12.M(),weight);
                  InvMass13Fake.at(ambiguity)->Fill(LV13.M(),weight);
                  InvMass23Fake.at(ambiguity)->Fill(LV23.M(),weight);
		}
	      }
	    }
	    const reco::PFCandidateRefVector & cands = thePFTau->signalPFChargedHadrCands();
	    if(cands.size()>=3 && thePFTau->decayMode()==10){
	      FakeRateHPS_eta.at(ambiguity)->Fill(thePFTau->eta(),weight);
              FakeRateHPS_pt.at(ambiguity)->Fill(thePFTau->pt(),weight);
	    }
	  }
	}
	         
	FakeRate_pt_isFit_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRate_pt_isFit_Eff.at(ambiguity)->getTH1F()->Divide(FakeRate_pt_isFit.at(ambiguity)->getTH1F(),FakeRate_pt_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	
	FakeRate_eta_isFit_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRate_eta_isFit_Eff.at(ambiguity)->getTH1F()->Divide(FakeRate_eta_isFit.at(ambiguity)->getTH1F(),FakeRate_eta_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	
	FakeRate_pt_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRate_pt_Eff.at(ambiguity)->getTH1F()->Divide(FakeRate_pt.at(ambiguity)->getTH1F(),FakeRate_pt_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	
	FakeRate_eta_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRate_eta_Eff.at(ambiguity)->getTH1F()->Divide(FakeRate_eta.at(ambiguity)->getTH1F(),FakeRate_eta_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	
	FakeRateHPS_pt_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRateHPS_pt_Eff.at(ambiguity)->getTH1F()->Divide(FakeRateHPS_pt.at(ambiguity)->getTH1F(),FakeRate_pt_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	
	FakeRateHPS_eta_Eff.at(ambiguity)->getTH1F()->Reset();
	FakeRateHPS_eta_Eff.at(ambiguity)->getTH1F()->Divide(FakeRateHPS_eta.at(ambiguity)->getTH1F(),FakeRate_eta_All.at(ambiguity)->getTH1F(), 1., 1., "b");
	}
      }
    }
  }
  if(doDQM_){
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(gensrc_, genParticles);

    /*
    edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > SecVertexInfo;
    iEvent.getByLabel("secondaryVertexTagInfos", SecVertexInfo);
    //std::cout << "N SecVertexInfo " << SecVertexInfo->size() << std::endl;
    for(unsigned int i=0;i<SecVertexInfo->size();i++){
      //std::cout << "SecVertexInfo " << i << " " << SecVertexInfo->at(i).nVertices() << std::endl;

      // & primaryVertex() 
      for(unsigned int j=0;j<SecVertexInfo->at(i).nVertices();j++){

	std::cout << "Found SecVertexInfo: " <<  SecVertexInfo->at(i).secondaryVertex(j).refittedTracks().size() << std::endl;
	//
	reco::Vertex V=SecVertexInfo->at(i).secondaryVertex(j);
        const std::vector<reco::Track> tracks=V.refittedTracks();
        double SumPx(0),SumPy(0),SumPz(0),SumE(0);
	int charge=0;
	TLorentzVector a1(0,0,0,0);
        for(unsigned int i=0; i<tracks.size(); i++){
          SumPx += tracks.at(i).px();
          SumPy += tracks.at(i).py();
          SumPz += tracks.at(i).pz();
          SumE  += sqrt(pow(tracks.at(i).p(),2)+pow(PDGInfo::pi_mass(),2));
	  charge+=tracks.at(i).charge();
	  TLorentzVector LV(tracks.at(i).px(),tracks.at(i).py(),tracks.at(i).pz(),sqrt(pow(tracks.at(i).p(),2)+pow(PDGInfo::pi_mass(),2)));
	  a1+=LV;
        }
        if(sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2))<PDGInfo::tau_mass() && abs(charge)==1){
          if(TMath::Prob(V.chi2(),V.ndof())>0.01){
	    if(genParticles.isValid()){
	      bool foundtruth=false;
	      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end() && !foundtruth; ++itr){
		const reco::GenParticle mytau=(*itr);
		if(Tau_JAKID_Filter::isTruthTauInAcceptance(mytau,TauPtMin_,TauEtaMax_)){
		  TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
		  TauDecay_CMSSWReco TD;
		  unsigned int jak_id, TauBitMask;
		  TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
		  std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
		  if(A1Matched(DecayProd,a1) && doJAKID(jak_id) &&  (int)TD.nProng(TauBitMask)==3){
		    TLorentzVector mc_a1(0,0,0,0);
		    for(unsigned int j=0;j<DecayProd.size();j++){
		      if(fabs(DecayProd.at(j)->pdgId())==PdtPdgMini::a_1_plus){
			mc_a1.SetPxPyPzE(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
		      }
		    }
		    for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
		      JAKIDHPS.at(ambiguity)->Fill(jak_id,weight);
		      a1MassHPS.at(ambiguity)->Fill(a1.M(),weight);
		      Truth_a1MassHPS.at(ambiguity).at(0)->Fill(a1.M()-mc_a1.M(),weight);
		    }
		  }
		}
	      }
	    }
	  }
        }
      }
    }
    */

    edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

    
    edm::Handle<reco::PFTauCollection> inputCollection;
    iEvent.getByLabel(tauType_+"PFTauProducer", inputCollection);
    for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
      reco::PFTauRef thePFTau(inputCollection, iPFTau);
      const reco::PFCandidateRefVector & cands = thePFTau->signalPFChargedHadrCands();  
      if(fabs(thePFTau->eta())<TauEtaMax_ && thePFTau->pt()>TauPtMin_){
	if(cands.size()>=3 && thePFTau->decayMode()==10){
	  TLorentzVector Pmax(1,0,0,0);
	  TLorentzVector a1(thePFTau->px(),thePFTau->py(),thePFTau->pz(),thePFTau->energy());
	  for(unsigned int j=0;j<cands.size();j++){
	    TLorentzVector LV;LV.SetXYZM(cands.at(j)->px(),cands.at(j)->py(),cands.at(j)->pz(),PDGInfo::pi_mass());
	    if(LV.Pt()>Pmax.Pt())Pmax=LV;
	  }
	  double drmax=0;
	  for(unsigned int j=0;j<cands.size();j++){
	    TLorentzVector LV;LV.SetXYZM(cands.at(j)->px(),cands.at(j)->py(),cands.at(j)->pz(),PDGInfo::pi_mass());
	    if(LV.DeltaR(Pmax)>drmax) drmax=LV.DeltaR(Pmax);
	  }
	
	  //std::cout << "N Signal: " << cands.size() << " " << drmax << std::endl;
	  if(genParticles.isValid()){
	    bool foundtruth=false;
	    for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end() && !foundtruth; ++itr){
	      const reco::GenParticle mytau=(*itr);
	      if(Tau_JAKID_Filter::isTruthTauInAcceptance(mytau,TauPtMin_,TauEtaMax_)){
		TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
		TauDecay_CMSSWReco TD;
		unsigned int jak_id, TauBitMask;
		TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
		std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
		if(A1Matched(DecayProd,a1) && doJAKID(jak_id) &&  (int)TD.nProng(TauBitMask)==3){
		  TLorentzVector mc_a1(0,0,0,0);
                  for(unsigned int j=0;j<DecayProd.size();j++){
                    if(fabs(DecayProd.at(j)->pdgId())==PdtPdgMini::a_1_plus){
                      mc_a1.SetPxPyPzE(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
                    }
                  }
		  for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
		    JAKIDHPS.at(ambiguity)->Fill(jak_id,weight);
		    a1MassHPS.at(ambiguity)->Fill(a1.M(),weight);
		    Truth_a1MassHPS.at(ambiguity).at(0)->Fill(a1.M()-mc_a1.M(),weight);

		  }
		}
	      }
	    }
	  }
	  for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){ 
	    PionDrHPS.at(ambiguity)->Fill(drmax,weight);

	    std::vector<reco::TrackRef> ITT;
	    const reco::PFCandidateRefVector &cands = thePFTau->signalPFChargedHadrCands();
	    for(unsigned int i=0;i<cands.size();i++){
	      if(cands.at(i)->trackRef().isNonnull()) {                                                                                                                                                     
		ITT.push_back(cands.at(i)->trackRef());         
	      }     
	    }   
	    double minchi2prob(1.0),TrkQuality(5),maxPt(0),minPt(999),HFrac(1.0);
	    for(unsigned int i=0;i<ITT.size();i++){
	      double trackprob=TMath::Prob(ITT.at(i)->chi2(),ITT.at(i)->ndof());
	      if(minchi2prob>trackprob)minchi2prob=trackprob;
	      if(HFrac>ITT.at(i)->validFraction())HFrac=ITT.at(i)->validFraction();
	      if(maxPt<ITT.at(i)->pt())maxPt=ITT.at(i)->pt();
	      if(minPt>ITT.at(i)->pt())minPt=ITT.at(i)->pt();
	      if(ITT.at(i)->quality(reco::TrackBase::confirmed)){       if(TrkQuality>3) TrkQuality=3;}
	      else if(ITT.at(i)->quality(reco::TrackBase::highPurity)){ if(TrkQuality>2) TrkQuality=2;}
	      else if(ITT.at(i)->quality(reco::TrackBase::tight)){      if(TrkQuality>1) TrkQuality=1;}
	      else if(ITT.at(i)->quality(reco::TrackBase::loose)){      if(TrkQuality>0) TrkQuality=0;}//loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4                                                                        
	      //if(ITT.quality(reco::TrackBase::goodIterative) && TrkQuality>4) TrkQuality=4;                                                                                                                                                
	    }

	    Trackchi2HPS.at(ambiguity)->Fill(minchi2prob,weight);
	    TrackQualityHPS.at(ambiguity)->Fill(TrkQuality,weight);
	    MaxTrackPtHPS.at(ambiguity)->Fill(maxPt,weight);
	    MinTrackPtHPS.at(ambiguity)->Fill(minPt,weight);
	    MinoverMaxTrackPtHPS.at(ambiguity)->Fill(minPt/maxPt,weight);
	    HitFracHPS.at(ambiguity)->Fill(HFrac,weight);

	    TransientVertex tmpVtx_;
	    std::vector<reco::TransientTrack> trks_;
	    for (std::vector<reco::TrackRef>::const_iterator iter=ITT.begin(); iter!=ITT.end(); ++iter){
	      trks_.push_back(transTrackBuilder_->build(**iter));
	    }
	    edm::Handle<reco::BeamSpot> beamSpot;
	    iEvent.getByLabel("offlineBeamSpot",beamSpot);
	    reco::BeamSpot thebeamSpot=*beamSpot;
	    if (!SecondaryVertexHelper::checkSecVtx(trks_,tmpVtx_,false,&thebeamSpot))continue;
	    reco::Vertex primaryVertexReFit=tmpVtx_;

	    Vertexchi2HPS.at(ambiguity)->Fill(TMath::Prob(primaryVertexReFit.chi2(),3),weight);

	  }
	}
      }
    }
    
  //bool found=false;
  cnt_++;
  for(SelectedKinematicDecayCollection::const_iterator kinFitTau=KinematicFitTaus->begin();kinFitTau!=KinematicFitTaus->end();kinFitTau++){
    for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
      unsigned int npassed=0;
      for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr){
	if(kinFitTau->discriminators(ambiguity).count(*discr)>0){
	  if(kinFitTau->discriminators(ambiguity).find(*discr)->second) npassed++;
	}
      }
      if(npassed==discriminators_.size()){
	cntFound_.at(ambiguity)+=1;
	//found=true;
	//std::cout << "KinematicTauAnalyzer::analyze A" << std::endl;
	SelectedKinematicDecay KFTau=(*kinFitTau);
	const TLorentzVector Tau=KFTau.Tau(ambiguity);
	const TMatrixDSym TauCov=KFTau.topParticle(ambiguity)->matrix();
	const TLorentzVector a1=KFTau.a1_p4(ambiguity);
	const TLorentzVector a1_initial=KFTau.Initial_a1_p4();
	const TLorentzVector Nu=KFTau.Neutrino(ambiguity);
	const std::vector<TLorentzVector> Pions=KFTau.Pions(ambiguity);

	const TLorentzVector Tau_initial=KFTau.InitialTauGuess(ambiguity);
	const TLorentzVector Nu_initial=KFTau.InitialNeutrinoGuess(ambiguity);
	const std::vector<TLorentzVector> Pions_initial=KFTau.InitialPions();
	
	reco::Vertex Pvtx=KFTau.InitialPrimaryVertexReFit();
	reco::Vertex Pvtx_initial=KFTau.InitialPrimaryVertexReFit();
	
	reco::Vertex Secvtx=KFTau.SecondaryVertex(ambiguity);
	reco::Vertex Secvtx_initial=KFTau.InitialSecondaryVertex();

	TVector3 FlightDir_initial=KFTau.InitialTauFlghtDirGuess(ambiguity);
	TVector3 FlightDir(Secvtx.position().x()-Pvtx.position().x(),
			   Secvtx.position().y()-Pvtx.position().y(),
			   Secvtx.position().z()-Pvtx.position().z());
	//std::cout << "KinematicTauAnalyzer::analyze B" << std::endl;
	nEvt.at(ambiguity)->Fill(0.5,weight);
	Truth_TauMatched.at(ambiguity)->Fill(0.0,weight);
	TauMass.at(ambiguity)->Fill(Tau.M(),weight);
	dTauMass.at(ambiguity)->Fill(Tau.M()-Tau_initial.M(),weight);
	TauPhi.at(ambiguity)->Fill(Tau.Phi(),weight);
	TauTheta.at(ambiguity)->Fill(Tau.Theta(),weight);
	TauE.at(ambiguity)->Fill(Tau.E(),weight);
	TauPt.at(ambiguity)->Fill(Tau.Pt(),weight);
	TauFlightDir.at(ambiguity)->Fill(FlightDir.Angle(Tau.Vect()),weight);
        TauFlightDirInitial.at(ambiguity)->Fill(FlightDir_initial.Angle(Tau_initial.Vect()),weight);

	//std::cout << "DQM " << Tau_initial.Px() << " " << Tau_initial.Py() << " " << Tau_initial.Pz() << " " << Tau_initial.E() 
	// << " Vector " << FlightDir_initial.X() << " " << FlightDir_initial.Y() << " " <<FlightDir_initial.Z() << std::endl;
	//std::cout << "KinematicTauAnalyzer::analyze C" << std::endl;
	GFAngle.at(ambiguity)->Fill(Tau.Angle(a1.Vect()),weight);
	TLorentzVector a1_TauRF=a1;
	a1_TauRF.Boost(-Tau.BoostVector());
	GFAngle_TauRF.at(ambiguity)->Fill(cos(Tau.Angle(a1_TauRF.Vect())),weight);

	TLorentzVector Pmax(1,0,0,0);
	for(unsigned int i=0; i<Pions.size();i++){
	  if(Pions.at(i).Pt()>Pmax.P())Pmax=Pions.at(i);
	  PionMass.at(ambiguity)->Fill(Pions.at(i).M(),weight);
	  dPionMass.at(ambiguity)->Fill(Pions.at(i).M()-Pions_initial.at(i).M(),weight);
	  PionPhiChange.at(ambiguity)->Fill(Pions.at(i).DeltaPhi(Pions_initial.at(i)),weight);
	  PionThetaChange.at(ambiguity)->Fill(fabs(Pions.at(i).Theta()-Pions_initial.at(i).Theta()),weight);
	  PionEChange.at(ambiguity)->Fill(Pions.at(i).E()-Pions_initial.at(i).E(),weight);
	}
	double dr=0;
	for(unsigned int i=0; i<Pions.size();i++){
	  if(fabs(Pions.at(i).Pt()-Pmax.P())>1.0){
	    if(dr<Pmax.DeltaR(Pions.at(i))) dr=Pmax.DeltaR(Pions.at(i));
	  }
	}
	PionDr.at(ambiguity)->Fill(dr,weight);
	const reco::PFCandidateRefVector &cands =  KFTau.PFTauRef()->signalPFChargedHadrCands();
	if(cands.size()>=3 && KFTau.PFTauRef()->decayMode()==10){
	  TLorentzVector Pmax(1,0,0,0);
	  for(unsigned int j=0;j<cands.size();j++){
	    TLorentzVector LV;LV.SetXYZM(cands.at(j)->px(),cands.at(j)->py(),cands.at(j)->pz(),PDGInfo::pi_mass());
	    if(LV.Pt()>Pmax.Pt())Pmax=LV;
	  }
	  double drmax=0;
	  for(unsigned int j=0;j<cands.size();j++){
	    TLorentzVector LV;LV.SetXYZM(cands.at(j)->px(),cands.at(j)->py(),cands.at(j)->pz(),PDGInfo::pi_mass());
	    if(LV.DeltaR(Pmax)>drmax) drmax=LV.DeltaR(Pmax);
	  }
	  //std::cout << "N Signal: " << cands.size() << " " << drmax << std::endl;
	  PionDrHPSwithFit.at(ambiguity)->Fill(drmax,weight);
	}

	if(Pions.size()==3){
	  std::vector<int> Charges=KFTau.Pions_Charge();
	  int charge=0;
	  for(unsigned int i=0;i<Pions.size();i++){
	    charge+=Charges.at(i);
	  }
	  TLorentzVector LV1, LV2, LV3;
	  bool has1=false;;
	  for(unsigned int i=0;i<Pions.size();i++){
	    if(has1==false && Charges.at(i)==charge){ LV1=Pions.at(i); has1=true;}
	    else if(Charges.at(i)==charge){ LV2=Pions.at(i);}
	    else{ LV3=Pions.at(i); }
	  }
	  TLorentzVector LV13=LV1+LV3;
	  TLorentzVector LV23=LV2+LV3;
	  TLorentzVector LV12=LV1+LV2;
	  if(LV13.M()>LV23.M()){
	    Dalitz.at(ambiguity)->Fill(LV13.M(),LV23.M(),weight);
	    InvMass12.at(ambiguity)->Fill(LV12.M(),weight);
	    InvMass13.at(ambiguity)->Fill(LV13.M(),weight);
	    InvMass23.at(ambiguity)->Fill(LV23.M(),weight);
	  }
	  else{
	    Dalitz.at(ambiguity)->Fill(LV23.M(),LV13.M(),weight);
	    InvMass12.at(ambiguity)->Fill(LV12.M(),weight);
	    InvMass13.at(ambiguity)->Fill(LV13.M(),weight);
	    InvMass23.at(ambiguity)->Fill(LV23.M(),weight);
	  }
	}

	//std::cout << "KinematicTauAnalyzer::analyze D" << std::endl;
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
        constraints.at(ambiguity)->Fill(KFTau.constraints(ambiguity),weight);
        csum.at(ambiguity)->Fill(KFTau.csum(ambiguity),weight);
	mincsum.at(ambiguity)->Fill(KFTau.mincsum(ambiguity),weight);

	chi2.at(ambiguity)->Fill(KFTau.chi2(ambiguity),weight);
        ndf.at(ambiguity)->Fill(KFTau.ndf(ambiguity),weight);
	chi2prob.at(ambiguity)->Fill(KFTau.chi2prob(ambiguity),weight);
	//std::cout << "chi2 " << KFTau.chi2(ambiguity) << " " << KFTau.ndf(ambiguity) << std::endl;
        chi2Vtx.at(ambiguity)->Fill(Pvtx.chi2(),weight);
        ndfVtx.at(ambiguity)->Fill(Pvtx.ndof(),weight);
        chi2probVtx.at(ambiguity)->Fill(TMath::Prob(Pvtx.chi2(),Pvtx.ndof()),weight);
	//std::cout << "PVT " << Pvtx.chi2() << " " << Pvtx.ndof() << " " << TMath::Prob(Pvtx.chi2(),Pvtx.ndof()) << std::endl;

	chi2SVtx.at(ambiguity)->Fill(KFTau.chi2Vtx(),weight);
	ndfSVtx.at(ambiguity)->Fill(KFTau.ndfVtx(),weight);
	chi2probSVtx.at(ambiguity)->Fill(TMath::Prob(KFTau.chi2Vtx(),KFTau.ndfVtx()),weight);
	//std::cout << "SVT " << KFTau.chi2Vtx() << " " << KFTau.ndfVtx() << " " << TMath::Prob(KFTau.chi2Vtx(),KFTau.ndfVtx()) << std::endl;

	FlightLength.at(ambiguity)->Fill(KFTau.FlightLength(),weight);
	FlightLengthSig.at(ambiguity)->Fill(KFTau.FlightLengthSig(),weight);	

	std::vector<reco::TrackRef> ITT=KFTau.InitialTrackTriplet();
	double minchi2prob(1.0),TrkQuality(5),maxPt(0),minPt(999),HFrac(1.0);
	for(unsigned int i=0;i<ITT.size();i++){
	  double trackprob=TMath::Prob(ITT.at(i)->chi2(),ITT.at(i)->ndof());
	  if(minchi2prob>trackprob)minchi2prob=trackprob;
	  if(HFrac>ITT.at(i)->validFraction())HFrac=ITT.at(i)->validFraction();
	  if(maxPt<ITT.at(i)->pt())maxPt=ITT.at(i)->pt();
	  if(minPt>ITT.at(i)->pt())minPt=ITT.at(i)->pt();
	  if(ITT.at(i)->quality(reco::TrackBase::confirmed)){       if(TrkQuality>3) TrkQuality=3;}
	  else if(ITT.at(i)->quality(reco::TrackBase::highPurity)){ if(TrkQuality>2) TrkQuality=2;}
	  else if(ITT.at(i)->quality(reco::TrackBase::tight)){      if(TrkQuality>1) TrkQuality=1;}
	  else if(ITT.at(i)->quality(reco::TrackBase::loose)){      if(TrkQuality>0) TrkQuality=0;}//loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4                                                              
	  //if(ITT.quality(reco::TrackBase::goodIterative) && TrkQuality>4) TrkQuality=4;                                                                                                                                      
	}

	Trackchi2.at(ambiguity)->Fill(minchi2prob,weight);
	TrackQuality.at(ambiguity)->Fill(TrkQuality,weight);
	MaxTrackPt.at(ambiguity)->Fill(maxPt,weight);
	MinTrackPt.at(ambiguity)->Fill(minPt,weight);
	MinoverMaxTrackPt.at(ambiguity)->Fill(minPt/maxPt,weight);
	HitFrac.at(ambiguity)->Fill(HFrac,weight);
	Vertexchi2.at(ambiguity)->Fill(TMath::Prob(KFTau.chi2Vtx(),3),weight);


	//std::cout << "KinematicTauAnalyzer::analyze E" << std::endl;
	// If Truth is valid run truth comparison
	if(genParticles.isValid()){
	  //std::cout << "KinematicTauAnalyzer::analyze F" << std::endl;
	  bool foundtruth=false;
	  for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end() && !foundtruth; ++itr){
	    const reco::GenParticle mytau=(*itr);
	    if(Tau_JAKID_Filter::isTruthTauInAcceptance(mytau,TauPtMin_,TauEtaMax_)){
	      //std::cout << "KinematicTauAnalyzer::analyze FA" << std::endl;
	      TLorentzVector mc(itr->p4().Px(),itr->p4().Py(),itr->p4().Pz(),itr->p4().E());
	      TauDecay_CMSSWReco TD;
	      unsigned int jak_id, TauBitMask;
	      TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	      std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
	      if(A1Matched(DecayProd,a1)){
		foundtruth=true;
		if(!(doJAKID(jak_id) &&  (int)TD.nProng(TauBitMask)==3)){
		  JAKIDFailed.at(ambiguity)->Fill(jak_id,weight);
		}
		else{
		  Truth_TauMatched.at(ambiguity)->Fill(1.0,weight);
		  //std::cout << "KinematicTauAnalyzer::analyze FC" << std::endl;
		  JAKID.at(ambiguity)->Fill(jak_id,weight);
		  //std::cout << "KinematicTauAnalyzer::analyze FD" << std::endl;
		  if(ambiguity==MultiProngTauSolver::plus || ambiguity==MultiProngTauSolver::minus){
		    TLorentzVector Tau_plus=KFTau.Tau(MultiProngTauSolver::plus);
		    TLorentzVector Tau_minus=KFTau.Tau(MultiProngTauSolver::minus);
		    if(fabs(mc.E()-Tau_plus.E())<fabs(mc.E()-Tau_minus.E()) && Tau_plus.E()!=Tau.E())  continue;
		    if(fabs(mc.E()-Tau_plus.E())>fabs(mc.E()-Tau_minus.E()) && Tau_minus.E()!=Tau.E()) continue;
		  }
		  //std::cout << "KinematicTauAnalyzer::analyze FE" << std::endl;
		  TVector3 TruthPvtx(itr->vx(),itr->vy(),itr->vz());
		  TVector3 TruthSvtx;
		  for( unsigned int dp=0;dp<DecayProd.size();dp++){
		    if(fabs(DecayProd.at(dp)->pdgId())!=fabs(PdtPdgMini::tau_minus))
		      TruthSvtx=TVector3(DecayProd.at(dp)->vx(),DecayProd.at(dp)->vy(),DecayProd.at(dp)->vz());
		  }

		  TVector3 dist=TruthSvtx-TruthPvtx;
		  double length=dist.Mag();
		  //std::cout << "KinematicTauAnalyzer::analyze F1" << std::endl;
		  
		  if(JAKIDtoIndex.count(jak_id)==1){
		    unsigned int idx=JAKIDtoIndex.find(jak_id)->second;
		    TMatrixDSym Pvtx_cov(3);
		    Pvtx_cov.ResizeTo(TMatrixDSym(3));
		    for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Pvtx_cov(i,j) = Pvtx.covariance(i,j);//diagonals are squares of sigmas
		    TMatrixDSym Svtx_cov(3);
		    //std::cout << "KinematicTauAnalyzer::analyze F2" << std::endl;
		    Svtx_cov.ResizeTo(TMatrixDSym(3));
		    for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Svtx_cov(i,j) = Secvtx.covariance(i,j);//diagonals are squares of sigmas
		    TVector3 Pvtx_point(Pvtx.position().x(),Pvtx.position().y(),Pvtx.position().z());
		    TVector3 Svtx_point(Secvtx.position().x(),Secvtx.position().y(),Secvtx.position().z());
		    //std::cout << "KinematicTauAnalyzer::analyze F3" << std::endl;
		    TMatrixDSym Truth_cov(3);
		    Truth_cov.ResizeTo(TMatrixDSym(3));
		    for(int i=0; i!=3; i++) for(int j=0; j!=3; j++) Svtx_cov(i,j) = Secvtx.covariance(i,j);
		    //std::cout << "KinematicTauAnalyzer::analyze G" << std::endl;
		    //Tau
		    Truth_TauMatch_dPhi.at(ambiguity).at(idx)->Fill(Tau.DeltaPhi(mc),weight);
		    Truth_TauMatch_dTheta.at(ambiguity).at(idx)->Fill(mc.Theta()-Tau.Theta(),weight);
		    Truth_TauMatch_dE.at(ambiguity).at(idx)->Fill(mc.E()-Tau.E(),weight);
		    Truth_TauMatch_dPt.at(ambiguity).at(idx)->Fill(mc.Pt()-Tau.Pt(),weight);
		    Truth_TauMatch_dPz.at(ambiguity).at(idx)->Fill(mc.Pz()-Tau.Pz(),weight);
		    
		    
		    TruthVtxX.at(ambiguity).at(idx)->Fill(Pvtx.position().x()-TruthPvtx.X(),weight);
		    TruthVtxY.at(ambiguity).at(idx)->Fill(Pvtx.position().y()-TruthPvtx.Y(),weight);
		    TruthVtxZ.at(ambiguity).at(idx)->Fill(Pvtx.position().z()-TruthPvtx.Z(),weight);
		  
                  PullVtxX.at(ambiguity).at(idx)->Fill((Pvtx.position().x()-TruthPvtx.X())/sqrt(Pvtx_cov(0,0)),weight);
		  PullVtxY.at(ambiguity).at(idx)->Fill((Pvtx.position().y()-TruthPvtx.Y())/sqrt(Pvtx_cov(1,1)),weight);
		  PullVtxZ.at(ambiguity).at(idx)->Fill((Pvtx.position().z()-TruthPvtx.Z())/sqrt(Pvtx_cov(2,2)),weight);

		  TruthSecVtxX.at(ambiguity).at(idx)->Fill(Secvtx.position().x()-TruthSvtx.X(),weight);
		  TruthSecVtxY.at(ambiguity).at(idx)->Fill(Secvtx.position().y()-TruthSvtx.Y(),weight);
		  TruthSecVtxZ.at(ambiguity).at(idx)->Fill(Secvtx.position().z()-TruthSvtx.Z(),weight);

                  PullSecVtxX.at(ambiguity).at(idx)->Fill((Secvtx.position().x()-TruthSvtx.X())/sqrt(Svtx_cov(0,0)),weight);
                  PullSecVtxY.at(ambiguity).at(idx)->Fill((Secvtx.position().y()-TruthSvtx.Y())/sqrt(Svtx_cov(1,1)),weight);
		  PullSecVtxZ.at(ambiguity).at(idx)->Fill((Secvtx.position().z()-TruthSvtx.Z())/sqrt(Svtx_cov(2,2)),weight);

		  /////////////////////////////////////
		  // rotate into tau truth direction
		  TMatrixT<double> Res(5,1);
		  Res(0,0)=Secvtx.position().x()-TruthSvtx.X();
		  Res(1,0)=Secvtx.position().y()-TruthSvtx.Y();
		  Res(2,0)=Secvtx.position().z()-TruthSvtx.Z();
		  Res(3,0)=mc.Phi();
		  Res(4,0)=mc.Theta();
		  TMatrixTSym<double> ResCov(5);
		  for(unsigned int s=0;s<3;s++){
		    for(unsigned int t=0;t<3;t++){
		      ResCov(s,t)=Svtx_cov(s,t);
		    }
		  }
		  TMatrixT<double> Resp=MultiProngTauSolver::RotateToTauFrame(Res);
		  TMatrixTSym<double> RespCov=ErrorMatrixPropagator::PropogateError(&MultiProngTauSolver::RotateToTauFrame,Res,ResCov);
                  TruthSecVtxXp.at(ambiguity).at(idx)->Fill(Resp(0,0),weight);
		  TruthSecVtxYp.at(ambiguity).at(idx)->Fill(Resp(1,0),weight);
                  TruthSecVtxZp.at(ambiguity).at(idx)->Fill(Resp(2,0),weight);

                  PullSecVtxXp.at(ambiguity).at(idx)->Fill(Resp(0,0)/sqrt(RespCov(0,0)),weight);
                  PullSecVtxYp.at(ambiguity).at(idx)->Fill(Resp(1,0)/sqrt(RespCov(1,1)),weight);
                  PullSecVtxZp.at(ambiguity).at(idx)->Fill(Resp(2,0)/sqrt(RespCov(2,2)),weight);

		  /////////////////////////////////////

		  if(TauCov(LorentzVectorParticle::px,LorentzVectorParticle::px)!=0 && TauCov(LorentzVectorParticle::py,LorentzVectorParticle::py)!=0 && TauCov(LorentzVectorParticle::pz,LorentzVectorParticle::pz)!=0){
		    PullTauPx.at(ambiguity).at(idx)->Fill((mc.Px()-Tau.Px())/sqrt(TauCov(LorentzVectorParticle::px,LorentzVectorParticle::px)),weight);
		    PullTauPy.at(ambiguity).at(idx)->Fill((mc.Py()-Tau.Py())/sqrt(TauCov(LorentzVectorParticle::py,LorentzVectorParticle::py)),weight);
		    PullTauPz.at(ambiguity).at(idx)->Fill((mc.Pz()-Tau.Pz())/sqrt(TauCov(LorentzVectorParticle::pz,LorentzVectorParticle::pz)),weight);
		  }

		  //std::cout << "Pull " << (Secvtx.position().x()-TruthSvtx.X())/sqrt(Svtx_cov(0,0)) << " Fit " << Secvtx.position().x() << " truth " << TruthSvtx.X() << " error " << sqrt(Svtx_cov(0,0)) << std::endl;
		  
		  VertexRotation vtxC;
		  TruthPVtxSig.at(ambiguity).at(idx)->Fill(vtxC.vtxDistanceSignificance(Pvtx_point,Pvtx_cov,TruthPvtx,Truth_cov),weight);
		  TruthSecVtxSig.at(ambiguity).at(idx)->Fill(vtxC.vtxDistanceSignificance(Svtx_point,Svtx_cov,TruthSvtx,Truth_cov),weight);

		  TVector3 TruthFlightDir=TruthSvtx-TruthPvtx;
		  TruthTauFlightDir.at(ambiguity).at(idx)->Fill(fabs(Tau.Angle(mc.Vect())),weight);
		  TruthTauFlightDirCheck.at(ambiguity).at(idx)->Fill(fabs(TruthFlightDir.Angle(mc.Vect())),weight);

                  Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.Pt()-mc.Pt(),weight);
		  Truth_TauMatch_reldPtvsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
                  Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.E()-mc.E(),weight);
		  Truth_TauMatch_dphivsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.DeltaPhi(mc),weight);
		  Truth_TauMatch_dthetavsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.Theta()-mc.Theta(),weight);

                  Truth_TauMatch_dPtvsTheta.at(ambiguity).at(idx)->Fill(mc.Theta(),Tau.Pt()-mc.Pt(),weight);
                  Truth_TauMatch_reldPtvsTheta.at(ambiguity).at(idx)->Fill(mc.Theta(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
                  Truth_TauMatch_dEvsTheta.at(ambiguity).at(idx)->Fill(mc.Theta(),Tau.E()-mc.E(),weight);
                  Truth_TauMatch_dphivsTheta.at(ambiguity).at(idx)->Fill(mc.Theta(),Tau.DeltaPhi(mc),weight);
                  Truth_TauMatch_dthetavsTheta.at(ambiguity).at(idx)->Fill(mc.Theta(),Tau.Theta()-mc.Theta(),weight);

                  Truth_TauMatch_dPtvsEta.at(ambiguity).at(idx)->Fill(mc.Eta(),Tau.Pt()-mc.Pt(),weight);
                  Truth_TauMatch_reldPtvsEta.at(ambiguity).at(idx)->Fill(mc.Eta(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
                  Truth_TauMatch_dEvsEta.at(ambiguity).at(idx)->Fill(mc.Eta(),Tau.E()-mc.E(),weight);
                  Truth_TauMatch_dphivsEta.at(ambiguity).at(idx)->Fill(mc.Eta(),Tau.DeltaPhi(mc),weight);
                  Truth_TauMatch_dthetavsEta.at(ambiguity).at(idx)->Fill(mc.Eta(),Tau.Theta()-mc.Theta(),weight);


		  Truth_TauMatch_dPtvsPt.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.Pt()-mc.Pt(),weight);
		  Truth_TauMatch_dEvsPt.at(ambiguity).at(idx)->Fill(mc.Pt(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
		  Truth_TauMatch_dphivsPt.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.E()-mc.E(),weight);
		  Truth_TauMatch_dthetavsPt.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.DeltaPhi(mc),weight);
		  Truth_TauMatch_reldPtvsPt.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.Theta()-mc.Theta(),weight);
		  if(fabs(KFTau.PFTauRef()->eta())<1.0){
		    Truth_TauMatch_dPtvsPtandEtaCut.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.Pt()-mc.Pt(),weight);
		    Truth_TauMatch_dEvsPtandEtaCut.at(ambiguity).at(idx)->Fill(mc.Pt(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
		    Truth_TauMatch_dphivsPtandEtaCut.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.E()-mc.E(),weight);
		    Truth_TauMatch_dthetavsPtandEtaCut.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.DeltaPhi(mc),weight);
		    Truth_TauMatch_reldPtvsPtandEtaCut.at(ambiguity).at(idx)->Fill(mc.Pt(),Tau.Theta()-mc.Theta(),weight);
		    
		    Truth_TauMatch_dPtvsLandEtaCut.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.Pt()-mc.Pt(),weight);
		    Truth_TauMatch_dEvsLandEtaCut.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
		    Truth_TauMatch_dphivsLandEtaCut.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.E()-mc.E(),weight);
		    Truth_TauMatch_dthetavsLandEtaCut.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.DeltaPhi(mc),weight);
		    Truth_TauMatch_reldPtvsLandEtaCut.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),Tau.Theta()-mc.Theta(),weight);
		    
		    Truth_TauMatch_dPtvsPhiandEtaCut.at(ambiguity).at(idx)->Fill(mc.Phi(),Tau.Pt()-mc.Pt(),weight);
		    Truth_TauMatch_dEvsPhiandEtaCut.at(ambiguity).at(idx)->Fill(mc.Phi(),(Tau.Pt()-mc.Pt())/mc.Pt(),weight);
		    Truth_TauMatch_dphivsPhiandEtaCut.at(ambiguity).at(idx)->Fill(mc.Phi(),Tau.E()-mc.E(),weight);
		    Truth_TauMatch_dthetavsPhiandEtaCut.at(ambiguity).at(idx)->Fill(mc.Phi(),Tau.DeltaPhi(mc),weight);
		    Truth_TauMatch_reldPtvsPhiandEtaCut.at(ambiguity).at(idx)->Fill(mc.Phi(),Tau.Theta()-mc.Theta(),weight);
		  }

		  //std::cout << "KinematicTauAnalyzer::analyze I" << std::endl;
		  TLorentzVector mc_a1(0,0,0,0);
		  for(unsigned int j=0;j<DecayProd.size();j++){
		    if(fabs(DecayProd.at(j)->pdgId())==PdtPdgMini::a_1_plus){
		      mc_a1.SetPxPyPzE(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
		    }
		  }
		  Truth_TauMatch_dGFAngle.at(ambiguity).at(idx)->Fill(a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->Fill(TruthFlightDir.Mag(),a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  
		  ///////////////////
		  //Quality cuts and GF analge
		  // track chi2
		  std::vector<reco::TrackRef> ITT=KFTau.InitialTrackTriplet();
		  double minchi2prob(1.0),TrkQuality(5),maxPt(0),minPt(999),HFrac(1.0);
		  for(unsigned int i=0;i<ITT.size();i++){
		    double trackprob=TMath::Prob(ITT.at(i)->chi2(),ITT.at(i)->ndof());
		    if(minchi2prob>trackprob)minchi2prob=trackprob;
		    if(HFrac>ITT.at(i)->validFraction())HFrac=ITT.at(i)->validFraction();
		    if(maxPt<ITT.at(i)->pt())maxPt=ITT.at(i)->pt();
		    if(minPt>ITT.at(i)->pt())minPt=ITT.at(i)->pt();
		    if(ITT.at(i)->quality(reco::TrackBase::confirmed)){       if(TrkQuality>3) TrkQuality=3;}
		    else if(ITT.at(i)->quality(reco::TrackBase::highPurity)){ if(TrkQuality>2) TrkQuality=2;}
		    else if(ITT.at(i)->quality(reco::TrackBase::tight)){      if(TrkQuality>1) TrkQuality=1;}
		    else if(ITT.at(i)->quality(reco::TrackBase::loose)){      if(TrkQuality>0) TrkQuality=0;}//loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4
		    //if(ITT.quality(reco::TrackBase::goodIterative) && TrkQuality>4) TrkQuality=4;
		  }
		  Truth_TauMatch_dGFAnglevsTrackchi2.at(ambiguity).at(idx)->Fill(minchi2prob,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
                  Truth_TauMatch_dGFAnglevsTrackQuality.at(ambiguity).at(idx)->Fill(TrkQuality,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFAnglevsMaxTrackPt.at(ambiguity).at(idx)->Fill(maxPt,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
                  Truth_TauMatch_dGFAnglevsMinTrackPt.at(ambiguity).at(idx)->Fill(minPt,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
                  Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt.at(ambiguity).at(idx)->Fill(minPt/maxPt,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFAnglevsHitFrac.at(ambiguity).at(idx)->Fill(HFrac,a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  Truth_TauMatch_dGFAnglevsVertexchi2.at(ambiguity).at(idx)->Fill(TMath::Prob(KFTau.chi2Vtx(),3),a1.Angle(Tau.Vect())-mc.Angle(mc_a1.Vect()),weight);
		  //std::cout << "Prob " << TMath::Prob(KFTau.chi2(ambiguity),3) << " " << KFTau.chi2(ambiguity) << std::endl;
		  ///////////////////

		  Truth_TauMatch_TrueMaxGFAngle.at(ambiguity).at(idx)->Fill(MultiProngTauSolver::ThetaGJMax(mc_a1),weight);
		  Truth_TauMatch_TrueGFAngle.at(ambiguity).at(idx)->Fill(mc.Angle(mc_a1.Vect()),weight);
		  TLorentzVector mc_a1_TauRF=mc_a1;
		  mc_a1_TauRF.Boost(-mc.BoostVector());
		  Truth_TauMatch_TrueGFAngle_TauRF.at(ambiguity).at(idx)->Fill(cos(mc.Angle(mc_a1_TauRF.Vect())),weight);
		  //
		  TLorentzVector a1_TauRF=a1;
		  a1_TauRF.Boost(-Tau.BoostVector());
		  Truth_TauMatch_GFAngleRecovsTruth_TauRF.at(ambiguity).at(idx)->Fill(cos(mc.Angle(mc_a1_TauRF.Vect())),cos(Tau.Angle(a1_TauRF.Vect())),weight);

		  //
		  if(MultiProngTauSolver::ThetaGJMax(mc_a1)!=0)Truth_TauMatch_TrueGFAngleoverMaxGFAngle.at(ambiguity).at(idx)->Fill(mc.Angle(mc_a1.Vect())/MultiProngTauSolver::ThetaGJMax(mc_a1),weight);
		  //std::cout << "Tau  " << ambiguity << " " << idx << std::endl;
		  Truth_A1Match_dPhi.at(ambiguity).at(idx)->Fill(a1.DeltaPhi(mc_a1));
		  Truth_A1Match_dTheta.at(ambiguity).at(idx)->Fill(a1.Theta()-mc_a1.Theta());
		  Truth_A1Match_dE.at(ambiguity).at(idx)->Fill(a1.E()-mc_a1.E());
		  Truth_A1Match_dPt.at(ambiguity).at(idx)->Fill(a1.Pt()-mc_a1.Pt());
		  Truth_A1Match_M.at(ambiguity).at(idx)->Fill(a1.M()-mc_a1.M());
		  //
		  Truth_A1Match_dPtvslength.at(ambiguity).at(idx)->Fill(length,a1.Pt()-mc_a1.Pt());
		  Truth_A1Match_dphivslength.at(ambiguity).at(idx)->Fill(length,a1.DeltaPhi(mc_a1));
		  Truth_A1Match_dthetavslength.at(ambiguity).at(idx)->Fill(length,a1.Theta()-mc_a1.Theta());
                  Truth_A1Match_reldPtvslength.at(ambiguity).at(idx)->Fill(length,(a1.Pt()-mc_a1.Pt())/mc_a1.Pt());

		  //std::cout << "KinematicTauAnalyzer::analyze J" << std::endl;
		  //charged hadrons (pi/K)
		  TLorentzVector Pmax(1,0,0,0);
		  for(unsigned int i=0; i<Pions.size();i++){
		    if(Pions.at(i).Pt()>Pmax.P())Pmax=Pions.at(i);
		  }

		  for(unsigned int i=0; i<Pions.size();i++){
		    double drmax=0;
		    double pidrmin=999;
		    TLorentzVector mcpion;
		    for(unsigned int j=0;j<DecayProd.size();j++){
		      if((fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::pi_plus) ||fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::K_plus)) && DecayProd.at(j)->status()==1){

			TLorentzVector mcPions_t(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
			double delta=sqrt(pow(mcPions_t.DeltaPhi(Pions.at(i)),2.0)+pow(mcPions_t.Theta()-Pions.at(i).Theta(),2.0));
			if(drmax<Pmax.DeltaR(mcPions_t)) drmax=Pmax.DeltaR(mcPions_t);
			if(delta<pidrmin){
			  mcpion=mcPions_t;
			  pidrmin=delta;
			}
		      }
		      Truth_MatchPionDr.at(ambiguity).at(idx)->Fill(pidrmin,weight);
		    }
		    if(i==0){
		      Truth_PionDr.at(ambiguity).at(idx)->Fill(drmax,weight);
		    }
		    if(pidrmin<0.1){
		      Truth_PionMatch_dPhi.at(ambiguity).at(idx)->Fill(mcpion.DeltaPhi(Pions.at(i)),weight);
		      Truth_PionMatch_dTheta.at(ambiguity).at(idx)->Fill(Pions.at(i).Theta()-mcpion.Theta(),weight);
		      Truth_PionMatch_dE.at(ambiguity).at(idx)->Fill(Pions.at(i).E()-mcpion.E(),weight);
		      //
		      Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->Fill(length,Pions.at(i).Pt()-mcpion.Pt());
		      Truth_PionMatch_dphivslength.at(ambiguity).at(idx)->Fill(length,Pions.at(i).DeltaPhi(mcpion));
		      Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->Fill(length,Pions.at(i).Theta()-mcpion.Theta());
		      Truth_PionMatch_reldPtvslength.at(ambiguity).at(idx)->Fill(length,(Pions.at(i).Pt()-mcpion.Pt())/mcpion.Pt());
		    }
		  }


		  //std::cout << "KinematicTauAnalyzer::analyze K" << std::endl;
		  //nu
		  TLorentzVector mcnu;
		  bool hasnu=false;
		  for(unsigned int i=0;i<DecayProd.size();i++){
		    if(fabs(PdtPdgMini::nu_tau)==fabs(DecayProd.at(i)->pdgId()) && DecayProd.at(i)->status()==1){ 
		      mcnu= TLorentzVector(DecayProd.at(i)->p4().Px(),DecayProd.at(i)->p4().Py(),DecayProd.at(i)->p4().Pz(),DecayProd.at(i)->p4().E()); hasnu=true;
		    }
		  }
		  if(hasnu){
		    Truth_NuMatch_dPhi.at(ambiguity).at(idx)->Fill(mcnu.DeltaPhi(Nu),weight);
		    Truth_NuMatch_dTheta.at(ambiguity).at(idx)->Fill(mcnu.Theta()-Nu.Theta(),weight);
		    Truth_NuMatch_dE.at(ambiguity).at(idx)->Fill(mcnu.E()-Nu.E(),weight);
		  }
		  }
		}
	      }
	    }
	  }
	  if(!foundtruth)JAKIDFailedMatch.at(ambiguity)->Fill(0.5,weight);
	}
      }
    }
  }
  // Fill truth info for all taus within acceptance

  if(genParticles.isValid()){
    for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
      bool hastau=false;
      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	const reco::GenParticle mytau=(*itr);
	if(Tau_JAKID_Filter::isTruthTauInAcceptance(mytau,TauPtMin_,TauEtaMax_)){
	  hastau=true;
	  TauDecay_CMSSWReco TD;
	  unsigned int jak_id(0), TauBitMask(0);
	  TD.AnalyzeTau(&mytau,jak_id,TauBitMask);
	  JAKIDall.at(ambiguity)->Fill(jak_id,weight);
	  //
	  if(jak_id==5){
	  std::vector<const reco::GenParticle* > DecayProd=TD.Get_TauDecayProducts();
	  std::vector<TLorentzVector> mcPions;
	  for(unsigned int j=0;j<DecayProd.size();j++){
	    if(fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::pi_plus) && DecayProd.at(j)->status()==1){
	      mcPions.push_back(TLorentzVector(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E()));
	    }
	  }
	  TLorentzVector Pmax(1,0,0,0);
	  for(unsigned int i=0;i<mcPions.size() && mcPions.size()>2;i++){
	    if(mcPions.at(i).Pt()>Pmax.Pt()) Pmax=mcPions.at(i);
	  }
	  double mydrmax=0;
	  TLorentzVector theTau(0,0,0,0);
	  for(unsigned int i=0;i<mcPions.size() && mcPions.size()>1;i++){
	    //std::cout << mcPions.at(i).Px() << " " << mcPions.at(i).Py() << " " << mcPions.at(i).Pz() << std::endl;
	    //std::cout << Pmax.Px() << " " << Pmax.Py() << " " << Pmax.Pz() << std::endl; 
	    if(mcPions.at(i).Pt()<Pmax.Pt()){
	      theTau+=mcPions.at(i);
	      //std::cout << mcPions.at(i).DeltaR(Pmax) << std::endl;
	      if(mcPions.at(i).DeltaR(Pmax)>mydrmax)mydrmax=mcPions.at(i).DeltaR(Pmax);
	    }
	  }
	  if( mcPions.size()>2)Truth_PionDrAll.at(ambiguity).at(0)->Fill(mydrmax,weight);
	  //std::cout << "Pmax " << Pmax.Pt() << " " << mydrmax << " " << mcPions.size() << " " << DecayProd.size() << " " << theTau.M() << std::endl; 
	  }
	}
      }
      if(hastau){
	// nasty hack for Eff - recompute eff every event with a tau
	JAKIDeff.at(ambiguity)->getTH1F()->Reset();
	JAKIDeff.at(ambiguity)->getTH1F()->Divide(JAKID.at(ambiguity)->getTH1F(),JAKIDall.at(ambiguity)->getTH1F(), 1., 1., "b");
	JAKIDHPSeff.at(ambiguity)->getTH1F()->Reset();
        JAKIDHPSeff.at(ambiguity)->getTH1F()->Divide(JAKIDHPS.at(ambiguity)->getTH1F(),JAKIDall.at(ambiguity)->getTH1F(), 1., 1., "b");
	JAKIDFailedEff.at(ambiguity)->getTH1F()->Reset();
	JAKIDFailedEff.at(ambiguity)->getTH1F()->Divide(JAKIDFailed.at(ambiguity)->getTH1F(),JAKIDall.at(ambiguity)->getTH1F(), 1., 1., "b");
	JAKIDFailedMatchEff.at(ambiguity)->getTH1F()->Reset();
	JAKIDFailedMatchEff.at(ambiguity)->getTH1F()->Divide(JAKIDFailedMatch.at(ambiguity)->getTH1F(),nEvt.at(ambiguity)->getTH1F(),1., 1., "b");
      }
    }
  }
    
  }
}



void KinematicTauAnalyzer::beginJob(){
  cnt_ = 0;
  cntFound_ = std::vector<int>(MultiProngTauSolver::NAmbiguity,0);

  if(dbe){
    for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
      TString amb="";
      if(ambiguity==MultiProngTauSolver::zero)  amb="ZeroAmbiguity";
      if(ambiguity==MultiProngTauSolver::plus)  amb="PlusSolution";
      if(ambiguity==MultiProngTauSolver::minus) amb="MinusSolution";
      // Number of analyzed events
      dbe->setCurrentFolder("KinematicFitTau/FakRate");
      FakeRateHPS_eta.push_back(dbe->book1D("FakeRateHPS_eta"+amb,"FakeRateHPS_eta"+amb, 50, -2.0, 2.0));  FakeRateHPS_eta.at(ambiguity)->setAxisTitle("#eta");
      FakeRateHPS_pt.push_back(dbe->book1D("FakeRateHPS_pt"+amb,"FakeRateHPS_pt"+amb, 50, 0.0, 100.0));  FakeRateHPS_pt.at(ambiguity)->setAxisTitle("P_{t}");

      FakeRate_eta.push_back(dbe->book1D("FakeRate_eta"+amb,"FakeRate_eta"+amb, 50, -2.0, 2.0));  FakeRate_eta.at(ambiguity)->setAxisTitle("#eta");
      FakeRate_pt.push_back(dbe->book1D("FakeRate_pt"+amb,"FakeRate_pt"+amb, 50, 0.0, 100.0));  FakeRate_pt.at(ambiguity)->setAxisTitle("P_{t}");


      FakeRate_eta_All.push_back(dbe->book1D("FakeRate_eta_All"+amb,"FakeRate_eta"+amb, 50, -2.0, 2.0));  FakeRate_eta_All.at(ambiguity)->setAxisTitle("#eta");
      FakeRate_pt_All.push_back(dbe->book1D("FakeRate_pt_All"+amb,"FakeRate_pt"+amb, 50, 0.0, 100.0));  FakeRate_pt_All.at(ambiguity)->setAxisTitle("P_{t}");

      FakeRate_eta_Eff.push_back(dbe->book1D("FakeRate_eta_Eff"+amb,"FakeRate_eta"+amb, 50, -2.0, 2.0));  FakeRate_eta_Eff.at(ambiguity)->setAxisTitle("#eta");
      FakeRate_pt_Eff.push_back(dbe->book1D("FakeRate_pt_Eff"+amb,"FakeRate_pt"+amb, 50, 0.0, 100.0));  FakeRate_pt_Eff.at(ambiguity)->setAxisTitle("P_{t}");

      FakeRateHPS_eta_Eff.push_back(dbe->book1D("FakeRateHPS_eta_Eff"+amb,"FakeRateHPS_eta"+amb, 50, -2.0, 2.0));  FakeRateHPS_eta_Eff.at(ambiguity)->setAxisTitle("#eta");
      FakeRateHPS_pt_Eff.push_back(dbe->book1D("FakeRateHPS_pt_Eff"+amb,"FakeRateHPS_pt"+amb, 50, 0.0, 100.0));  FakeRateHPS_pt_Eff.at(ambiguity)->setAxisTitle("P_{t}");

      FakeRate_eta_isFit.push_back(dbe->book1D("FakeRate_eta_isFit"+amb,"FakeRate_eta"+amb, 50, -2.0, 2.0));  FakeRate_eta_isFit.at(ambiguity)->setAxisTitle("#eta");
      FakeRate_pt_isFit.push_back(dbe->book1D("FakeRate_pt_isFit"+amb,"FakeRate_pt"+amb, 50, 0.0, 100.0));  FakeRate_pt_isFit.at(ambiguity)->setAxisTitle("P_{t}");

      FakeRate_eta_isFit_Eff.push_back(dbe->book1D("FakeRate_eta_isFit_Eff"+amb,"FakeRate_eta"+amb, 50, -2.0, 2.0));  FakeRate_eta_isFit_Eff.at(ambiguity)->setAxisTitle("#eta");
      FakeRate_pt_isFit_Eff.push_back(dbe->book1D("FakeRate_pt_isFit_Eff"+amb,"FakeRate_pt"+amb, 50, 0.0, 100.0));  FakeRate_pt_isFit_Eff.at(ambiguity)->setAxisTitle("P_{t}");

      DalitzFake.push_back(dbe->book2D("DalitzFake"+amb,"DalitzFake",100,0.0,2.0,100,0.0,2.0));
      DalitzFake.at(ambiguity)->setAxisTitle("M_{#pi^{+}#pi^{-}} (GeV)",1); DalitzFake.at(ambiguity)->setAxisTitle("M_{#pi^{+}#pi^{-}} (GeV)",2);

      InvMass12Fake.push_back(dbe->book1D("InvMass12"+amb,"InvMass12",100,0.0,2.0));  InvMass12Fake.at(ambiguity)->setAxisTitle("M_{12} (GeV)");
      InvMass13Fake.push_back(dbe->book1D("InvMass13"+amb,"InvMass13",100,0.0,2.0));  InvMass13Fake.at(ambiguity)->setAxisTitle("M_{13} (GeV)");
      InvMass23Fake.push_back(dbe->book1D("InvMass23"+amb,"InvMass23",100,0.0,2.0));  InvMass23Fake.at(ambiguity)->setAxisTitle("M_{23} (GeV)");

      dbe->setCurrentFolder("KinematicFitTau/FitResult");
      nEvt.push_back(dbe->book1D("nEvt"+amb,"n analyzed Events "+amb, 1, 0., 1.));  nEvt.at(ambiguity)->setAxisTitle("Number of Events");
      TauMass.push_back(dbe->book1D("TauMass"+amb,"M_{Tau} "+amb,100,1.7,1.8));      TauMass.at(ambiguity)->setAxisTitle("M_{Tau} (GeV)");
      PionMass.push_back(dbe->book1D("PionMass"+amb,"M_{#pi} "+amb,100,0.13,0.14));  PionMass.at(ambiguity)->setAxisTitle("M_{#pi} (GeV)");
      NuMass.push_back(dbe->book1D("NuMass"+amb,"M_{nu} "+amb,100,-0.05,0.05));      NuMass.at(ambiguity)->setAxisTitle("M_{#nu} (GeV)");
      TauFlightDir.push_back(dbe->book1D("TauFlightDir"+amb,"|#psi_{#tau Dir.,Vtx}^{KF}-#psi_{#tau Dir.,#tau}^{KF}|"+amb,100,0.0,0.5));  TauFlightDir.at(ambiguity)->setAxisTitle("|#psi_{Vtx}^{#tau Dir.,KF}-#psi_{#tau Dir.,#tau}^{KF}| (rad)");
      TauFlightDirInitial.push_back(dbe->book1D("TauFlightDirInitial"+amb,"|#psi_{Vtx}^{Initial}-#psi_{#tau}^{Initial}|"+amb,100,0.0,0.5));  TauFlightDirInitial.at(ambiguity)->setAxisTitle("|#psi_{#tau Dir.,Vtx}^{Initial}-#psi_{#tau Dir.,#tau}^{Initial}| (rad)");

      ////////////////////////////////////////////////////////////////////
      dbe->setCurrentFolder("KinematicFitTau/InitialtoFit");
      VtxXChange.push_back(dbe->book1D("VtxXChange"+amb,"Vtx_{X} "+amb,100,-0.5,0.5));       VtxXChange.at(ambiguity)->setAxisTitle("Vtx_{X} (cm)"); 
      VtxYChange.push_back(dbe->book1D("VtxYChange"+amb,"Vtx_{Y} "+amb,100,-0.5,0.5));       VtxYChange.at(ambiguity)->setAxisTitle("Vtx_{Y} (cm)");
      VtxZChange.push_back(dbe->book1D("VtxZChange"+amb,"Vtx_{Z} "+amb,100,-0.5,0.5));       VtxZChange.at(ambiguity)->setAxisTitle("Vtx_{Z} (cm)");
      SecVtxXChange.push_back(dbe->book1D("SecVtxXChange"+amb,"Vtx_{X}^{Sec} "+amb,100,-0.02,0.02)); SecVtxXChange.at(ambiguity)->setAxisTitle("Vtx_{X}^{Sec} (cm)");
      SecVtxYChange.push_back(dbe->book1D("SecVtxYChange"+amb,"Vtx_{Y}^{Sec} "+amb,100,-0.02,0.02)); SecVtxYChange.at(ambiguity)->setAxisTitle("Vtx_{Y}^{Sec} (cm)");
      SecVtxZChange.push_back(dbe->book1D("SecVtxZChange"+amb,"Vtx_{Z}^{Sec} "+amb,100,-0.02,0.02)); SecVtxZChange.at(ambiguity)->setAxisTitle("Vtx_{Z}^{Sec} (cm)");
      
      TauPhi.push_back(dbe->book1D("TauPhi"+amb,"phi_{#tau} "+amb,100,-TMath::Pi(),TMath::Pi()));         TauPhi.at(ambiguity)->setAxisTitle("#phi_{#tau} (rad)");
      TauTheta.push_back(dbe->book1D("TauTheta"+amb,"theta_{#tau} "+amb,100,0.0,TMath::Pi()));  TauTheta.at(ambiguity)->setAxisTitle("#theta_{#tau} (rad)");
      TauE.push_back(dbe->book1D("TauE"+amb,"E_{#tau} "+amb,100,0,100));                   TauE.at(ambiguity)->setAxisTitle("E_{#tau} (GeV)");
      TauPt.push_back(dbe->book1D("TauPt"+amb,"P_{t}^{#tau} "+amb,100,0,100)); TauPt.at(ambiguity)->setAxisTitle("P_{t}^{#tau} (GeV)");
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
      a1Mass.push_back(dbe->book1D("a1Mass"+amb,"M_{a1}"+amb,100,0.0,2.5));  a1Mass.at(ambiguity)->setAxisTitle("M_{a1} (GeV)");
      energyTFraction.push_back(dbe->book1D("energyTFraction"+amb,"energyTFraction"+amb,100,0.0,2.0));  energyTFraction.at(ambiguity)->setAxisTitle("E_{a1}/E_{#tau}");
      iterations.push_back(dbe->book1D("iterations"+amb,"iterations"+amb,51,-0.5,50.5));  iterations.at(ambiguity)->setAxisTitle("Number of Iterations");
      maxiterations.push_back(dbe->book1D("maxiterations"+amb,"maxiterations"+amb,51,0.5,50.5));  maxiterations.at(ambiguity)->setAxisTitle("Max Number of Iterations");
      chi2.push_back(dbe->book1D("chi2"+amb,"LC chi2"+amb,100,0.0,25.0));  chi2.at(ambiguity)->setAxisTitle("LC#chi^{2}");
      constraints.push_back(dbe->book1D("constraints"+amb,"constraints"+amb,51,-0.5,50.5));  constraints.at(ambiguity)->setAxisTitle("Number of Constraints");
      ndf.push_back(dbe->book1D("ndf"+amb,"LC ndf"+amb,51,-0.5,50.5));  ndf.at(ambiguity)->setAxisTitle("LC N.D.F");
      csum.push_back(dbe->book1D("csum"+amb,"csum"+amb,51,-0.5,50.5));  csum.at(ambiguity)->setAxisTitle("csum");
      mincsum.push_back(dbe->book1D("mincsum"+amb,"mincsum"+amb,51,-0.5,50.5));  mincsum.at(ambiguity)->setAxisTitle("mincsum");
      chi2prob.push_back(dbe->book1D("chi2prob"+amb,"LC chi2prob"+amb,100,0.0,1.0));  chi2prob.at(ambiguity)->setAxisTitle("LC #chi^{2} Probabilty");
      GFAngle_TauRF.push_back(dbe->book1D("GFAngle_TauRF"+amb,"cos(#theta_{GF}^{#tau RF})"+amb,100,-1.0,1.0));  GFAngle_TauRF.at(ambiguity)->setAxisTitle("cos(#theta_{GF}^{#tau RF})");
      GFAngle.push_back(dbe->book1D("GFAngleKF"+amb,"#theta_{GF}^{Lab,KF}"+amb,100,0.0,0.5));  GFAngle.at(ambiguity)->setAxisTitle("#theta_{GF}^{KF,Initial} (rad)");


      chi2SVtx.push_back(dbe->book1D("chi2SVtx"+amb,"Secondary Vertex chi2"+amb,100,0.0,10.0));  chi2SVtx.at(ambiguity)->setAxisTitle("Secondary Vertex #chi^{2}");
      ndfSVtx.push_back(dbe->book1D("ndfSVtx"+amb,"Secondary Vertex ndf"+amb,51,-0.5,50.5));    ndfSVtx.at(ambiguity)->setAxisTitle("Secondary Vertex N.D.F");
      chi2probSVtx.push_back(dbe->book1D("chi2probSVtx"+amb,"Secondary Vertex chi2prob"+amb,100,0.0,1.0));  chi2probSVtx.at(ambiguity)->setAxisTitle("Secondary Vertex #chi^{2} Probabilty");
      chi2Vtx.push_back(dbe->book1D("chi2PVtx"+amb,"Primary Vertex chi2"+amb,100,0.0,10.0));  chi2Vtx.at(ambiguity)->setAxisTitle("Primary Vertex #chi^{2}");
      ndfVtx.push_back(dbe->book1D("ndfPVtx"+amb,"Primary Vertex ndf"+amb,51,-0.5,50.5));    ndfVtx.at(ambiguity)->setAxisTitle("Primary Vertex N.D.F");
      chi2probVtx.push_back(dbe->book1D("chi2probPVtx"+amb,"Primary Vertex chi2prob"+amb,100,0.0,1.0));  chi2probVtx.at(ambiguity)->setAxisTitle("Primary Vertex #chi^{2} Probabilty");
      FlightLength.push_back(dbe->book1D("FlightLength"+amb,"Flight-Length "+amb,100,0.0,5.0));  FlightLength.at(ambiguity)->setAxisTitle("Flight-Length (cm)");
      FlightLengthSig.push_back(dbe->book1D("FlightLengthSig"+amb,"Flight-Length Significance "+amb,200,0.0,50.0));  FlightLengthSig.at(ambiguity)->setAxisTitle("Flight-Length Significance");


      PionDr.push_back(dbe->book1D("PionDr"+amb,"Pion Dr",100,0.0,0.5));  PionDr.at(ambiguity)->setAxisTitle("Pion Dr_{max}");
      PionDrHPS.push_back(dbe->book1D("PionDrHPS"+amb,"Pion Dr",100,0.0,0.5));  PionDrHPS.at(ambiguity)->setAxisTitle("Pion Dr_{max}^{HPS}");
      PionDrHPSwithFit.push_back(dbe->book1D("PionDrHPSwithFit"+amb,"Pion Dr Sel",100,0.0,0.5));  PionDrHPSwithFit.at(ambiguity)->setAxisTitle("Pion Dr_{max}^{HPS,Sel}");

      a1MassHPS.push_back(dbe->book1D("a1MassHPS"+amb,"a1MassHPS",100,0.0,2.5));  a1MassHPS.at(ambiguity)->setAxisTitle("M_{a1}^{HPS}");
      
      Dalitz.push_back(dbe->book2D("Dalitz"+amb,"Dalitz",100,0.0,2.0,100,0.0,2.0));  Dalitz.at(ambiguity)->setAxisTitle("M_{#pi^{+}#pi^{-}} (GeV)",1); Dalitz.at(ambiguity)->setAxisTitle("M_{#pi^{+}#pi^{-}} (GeV)",2);  

      InvMass12.push_back(dbe->book1D("InvMass12"+amb,"InvMass12",100,0.0,2.0));  InvMass12.at(ambiguity)->setAxisTitle("M_{12} (GeV)");
      InvMass13.push_back(dbe->book1D("InvMass13"+amb,"InvMass13",100,0.0,2.0));  InvMass13.at(ambiguity)->setAxisTitle("M_{13} (GeV)");
      InvMass23.push_back(dbe->book1D("InvMass23"+amb,"InvMass23",100,0.0,2.0));  InvMass23.at(ambiguity)->setAxisTitle("M_{23} (GeV)");

      Trackchi2.push_back(dbe->book1D("Trackchi2"+amb,"prob( "+amb+")",20,0.0,1.0)); Trackchi2.at(ambiguity)->setAxisTitle("Track Probability",1);
      TrackQuality.push_back(dbe->book1D("TrackQuality"+amb,"quality( "+amb+")",5,-0.5,4.5)); TrackQuality.at(ambiguity)->setAxisTitle("Track Quaility",1);
      MaxTrackPt.push_back(dbe->book1D("MaxTrackPt"+amb,"P_{t}^{max}( "+amb+")",100,0.0,25.0)); MaxTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{max} (GeV)",1);
      MinTrackPt.push_back(dbe->book1D("MinTrackPt"+amb,"P_{t}^{min}( "+amb+")",40,0.0,10.0)); MinTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{min} (GeV)",1);
      MinoverMaxTrackPt.push_back(dbe->book1D("MinoverMaxTrackPt"+amb," P_{t}^{min}/P_{t}^{max}( "+amb+")",20,0.0,1.0));  MinoverMaxTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{min}/P_{t}^{max}",1);
      HitFrac.push_back(dbe->book1D("HitFrac"+amb,"Hit fraction( "+amb+")",22,0.0,1.1)); HitFrac.at(ambiguity)->setAxisTitle("HitFrac",1);
      Vertexchi2.push_back(dbe->book1D("Vertexchi2"+amb,"prob( "+amb+")",20,0.0,1.0)); Vertexchi2.at(ambiguity)->setAxisTitle("Vertex Probability",1);

      Trackchi2HPS.push_back(dbe->book1D("Trackchi2HPS"+amb,"prob( "+amb+")",20,0.0,1.0)); Trackchi2.at(ambiguity)->setAxisTitle("Track Probability",1);
      TrackQualityHPS.push_back(dbe->book1D("TrackQualityHPS"+amb,"quality( "+amb+")",5,-0.5,4.5)); TrackQuality.at(ambiguity)->setAxisTitle("Track Quaility",1);
      MaxTrackPtHPS.push_back(dbe->book1D("MaxTrackPtHPS"+amb,"P_{t}^{max}( "+amb+")",100,0.0,25.0)); MaxTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{max} (GeV)",1);
      MinTrackPtHPS.push_back(dbe->book1D("MinTrackPtHPS"+amb,"P_{t}^{min}( "+amb+")",40,0.0,10.0)); MinTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{min} (GeV)",1);
      MinoverMaxTrackPtHPS.push_back(dbe->book1D("MinoverMaxTrackPtHPS"+amb," P_{t}^{min}/P_{t}^{max}( "+amb+")",20,0.0,1.0));  MinoverMaxTrackPt.at(ambiguity)->setAxisTitle("P_{t}^{min}/P_{t}^{max}",1);
      HitFracHPS.push_back(dbe->book1D("HitFracHPS"+amb,"Hit fraction( "+amb+")",22,0.0,1.1)); HitFrac.at(ambiguity)->setAxisTitle("HitFrac",1);
      Vertexchi2HPS.push_back(dbe->book1D("Vertexchi2HPS"+amb,"prob( "+amb+")",20,0.0,1.0)); Vertexchi2.at(ambiguity)->setAxisTitle("Vertex Probability",1);


      /////////////////////////////////////////////////////////////////////
      // now for Truth
      dbe->setCurrentFolder("KinematicFitTau/Truth");      
      JAKID.push_back(dbe->book1D("JAKID"+amb,"JAK ID "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));              JAKID.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDHPS.push_back(dbe->book1D("JAKIDHPS"+amb,"JAK ID "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));        JAKIDHPS.at(ambiguity)->setAxisTitle("JAK ID HPS");
      JAKIDHPSeff.push_back(dbe->book1D("JAKIDHPSeff"+amb,"JAK ID "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));  JAKIDHPSeff.at(ambiguity)->setAxisTitle("JAK ID HPS Eff");
      JAKIDall.push_back(dbe->book1D("JAKIDall"+amb,"JAK ID All "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));    JAKIDall.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDeff.push_back(dbe->book1D("JAKIDeff"+amb,"JAK ID Eff "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));    JAKIDeff.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDFailed.push_back(dbe->book1D("JAKIDFailed"+amb,"JAK ID Eff "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));  JAKIDFailed.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDFailedEff.push_back(dbe->book1D("JAKIDFailedEff"+amb,"JAK ID Eff "+amb,TauDecay::NJAKID,-0.5,(float)(TauDecay::NJAKID)-0.5));  JAKIDFailedEff.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDFailedMatch.push_back(dbe->book1D("JAKIDFailedMatch"+amb,"JAK ID Eff "+amb,1, 0., 1.));        JAKIDFailedMatch.at(ambiguity)->setAxisTitle("JAK ID");
      JAKIDFailedMatchEff.push_back(dbe->book1D("JAKIDFailedMatchEff"+amb,"JAK ID Eff "+amb,1, 0., 1.));  JAKIDFailedMatchEff.at(ambiguity)->setAxisTitle("JAK ID");


      Truth_TauMatched.push_back(dbe->book1D("Truth_TauMatched"+amb,"Truth_TauMatched "+amb,2,-0.5,1.5));                       Truth_TauMatched.at(ambiguity)->setAxisTitle("Tau Matched (0=All 1=Matched)");



      Truth_TauMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dE.push_back(std::vector<MonitorElement*>());
 
      Truth_A1Match_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_dE.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_dPt.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_M.push_back(std::vector<MonitorElement*>());

      Truth_A1Match_dPtvslength.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_dphivslength.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_dthetavslength.push_back(std::vector<MonitorElement*>());
      Truth_A1Match_reldPtvslength.push_back(std::vector<MonitorElement*>());

      Truth_PionMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dE.push_back(std::vector<MonitorElement*>());

      Truth_PionMatch_dPtvslength.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dphivslength.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_dthetavslength.push_back(std::vector<MonitorElement*>());
      Truth_PionMatch_reldPtvslength.push_back(std::vector<MonitorElement*>());

      Truth_NuMatch_dPhi.push_back(std::vector<MonitorElement*>());
      Truth_NuMatch_dTheta.push_back(std::vector<MonitorElement*>());
      Truth_NuMatch_dE.push_back(std::vector<MonitorElement*>());

      TruthVtxX.push_back(std::vector<MonitorElement*>());
      TruthVtxY.push_back(std::vector<MonitorElement*>());
      TruthVtxZ.push_back(std::vector<MonitorElement*>());
      TruthSecVtxX.push_back(std::vector<MonitorElement*>());
      TruthSecVtxY.push_back(std::vector<MonitorElement*>());
      TruthSecVtxZ.push_back(std::vector<MonitorElement*>());
      TruthSecVtxXp.push_back(std::vector<MonitorElement*>());
      TruthSecVtxYp.push_back(std::vector<MonitorElement*>());
      TruthSecVtxZp.push_back(std::vector<MonitorElement*>());

      PullVtxX.push_back(std::vector<MonitorElement*>());
      PullVtxY.push_back(std::vector<MonitorElement*>());
      PullVtxZ.push_back(std::vector<MonitorElement*>());
      PullSecVtxX.push_back(std::vector<MonitorElement*>());
      PullSecVtxY.push_back(std::vector<MonitorElement*>());
      PullSecVtxZ.push_back(std::vector<MonitorElement*>());
      PullSecVtxXp.push_back(std::vector<MonitorElement*>());
      PullSecVtxYp.push_back(std::vector<MonitorElement*>());
      PullSecVtxZp.push_back(std::vector<MonitorElement*>());

      PullTauPx.push_back(std::vector<MonitorElement*>());
      PullTauPy.push_back(std::vector<MonitorElement*>());
      PullTauPz.push_back(std::vector<MonitorElement*>());

      TruthPVtxSig.push_back(std::vector<MonitorElement*>());
      TruthSecVtxSig.push_back(std::vector<MonitorElement*>());
      TruthTauFlightDir.push_back(std::vector<MonitorElement*>());
      TruthTauFlightDirCheck.push_back(std::vector<MonitorElement*>());

      Truth_TauMatch_dPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPz.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPtvsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsL.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAngle.push_back(std::vector<MonitorElement*>());

      Truth_TauMatch_dPtvsTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsTheta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsTheta.push_back(std::vector<MonitorElement*>());

      Truth_TauMatch_dPtvsEta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsEta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsEta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsEta.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsEta.push_back(std::vector<MonitorElement*>());

      //
      Truth_TauMatch_dPtvsPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPtvsPtandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsPtandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsPtandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsPtandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsPtandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPtvsLandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsLandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsLandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsLandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsLandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dPtvsPhiandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dEvsPhiandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dphivsPhiandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dthetavsPhiandEtaCut.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_reldPtvsPhiandEtaCut.push_back(std::vector<MonitorElement*>());
      //

      Truth_TauMatch_TrueMaxGFAngle.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_TrueGFAngle.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_TrueGFAngle_TauRF.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_GFAngleRecovsTruth_TauRF.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_TrueGFAngleoverMaxGFAngle.push_back(std::vector<MonitorElement*>());

      Truth_TauMatch_dGFAnglevsTrackchi2.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsTrackQuality.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsMaxTrackPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsMinTrackPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsHitFrac.push_back(std::vector<MonitorElement*>());
      Truth_TauMatch_dGFAnglevsVertexchi2.push_back(std::vector<MonitorElement*>());
      Truth_PionDr.push_back(std::vector<MonitorElement*>());
      Truth_PionDrAll.push_back(std::vector<MonitorElement*>());
      Truth_a1MassHPS.push_back(std::vector<MonitorElement*>());
      Truth_MatchPionDr.push_back(std::vector<MonitorElement*>());
      Truth_MatchPionDrHPS.push_back(std::vector<MonitorElement*>());

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
	  Truth_NuMatch_dPhi.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdPhi"+tmp+amb,"d#phi_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.40,0.40));           axis="d#phi_{"+tmp+"}^{#nu Match} (rad)"; Truth_NuMatch_dPhi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_NuMatch_dTheta.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdTheta"+tmp+amb, "d#theta_{"+tmp+"}^{#nu Match} "+amb,100 ,-0.40,0.40));    axis="d#theta_{"+tmp+"}^{#nu Match} (rad)"; Truth_NuMatch_dTheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_NuMatch_dE.at(ambiguity).push_back(dbe->book1D("TruthNuMatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{#nu Match} "+amb,100 ,-50.0,50.0));                      axis="dE_{"+tmp+"}^{#nu Match} (GeV)"; Truth_NuMatch_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());



	  Truth_A1Match_dPhi.at(ambiguity).push_back(dbe->book1D("TruthA1MatchdPhi"+tmp+amb,"d#phi_{"+tmp+"}^{a_{1} Match} "+amb,100 ,-0.02,0.02));       axis="d#phi_{"+tmp+"}^{a_{1} Match} (rad)"; Truth_A1Match_dPhi.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_A1Match_dTheta.at(ambiguity).push_back(dbe->book1D("TruthA1MatchdTheta"+tmp+amb,"d#theta_{"+tmp+"}^{a_{1} Match} "+amb,100 ,-0.04,0.04)); axis="d#theta_{"+tmp+"}^{a_{1} Match} (rad)"; Truth_A1Match_dTheta.at(ambiguity).at(idx)->setAxisTitle(axis.Data()); 
	  Truth_A1Match_dE.at(ambiguity).push_back(dbe->book1D("TruthA1MatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{a_{1} Match} "+amb,100 ,-15.0,15.0));       axis="dE_{"+tmp+"}^{a_{1} Match} (GeV)"; Truth_A1Match_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_A1Match_dPt.at(ambiguity).push_back(dbe->book1D("TruthA1MatchdEnergy"+tmp+amb,"dE_{"+tmp+"}^{a_{1} Match} "+amb,100 ,-15.0,15.0));       axis="dE_{"+tmp+"}^{a_{1} Match} (GeV)"; Truth_A1Match_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_A1Match_M.at(ambiguity).push_back(dbe->book1D("TruthA1MatchdM"+tmp+amb,"dM_{"+tmp+"}^{a_{1} Match} "+amb,200 ,-0.5,0.5));       axis="dM_{"+tmp+"}^{a_{1} Match} (GeV)"; Truth_A1Match_dE.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	 
	  Truth_A1Match_dPtvslength.at(ambiguity).push_back(dbe->book2D("Truth_A1Match_dPtvslength"+tmp+amb,"Truth A1 Tau-Matchd dP_{t}^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,100,-5,5)); Truth_A1Match_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("P_{T}-P_{T}^{Truth} (GeV)",2); Truth_A1Match_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_A1Match_dphivslength.at(ambiguity).push_back(dbe->book2D("Truth_A1Match_dphivslength"+tmp+amb,"Truth A1 Tau-Matchd d#phi^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,100,-0.02,0.02)); Truth_A1Match_dphivslength.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (rad)",2); Truth_A1Match_dphivslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_A1Match_dthetavslength.at(ambiguity).push_back(dbe->book2D("Truth_A1Match_dthetavslength"+tmp+amb,"Truth A1 Tau-Matchd d#theta^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,100,-0.02,0.02)); Truth_A1Match_dthetavslength.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{Truth} (rad)",2); Truth_A1Match_dthetavslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_A1Match_reldPtvslength.at(ambiguity).push_back(dbe->book2D("Truth_A1Match_reldPtvslength"+tmp+amb,"Truth A1 Tau-Matchd Relative dP_{t}^{Lab}vs L( "+tmp+amb+")",100,0.0,5.0,100,-5,5)); Truth_A1Match_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_A1Match_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);

	  Truth_PionMatch_dPtvslength.at(ambiguity).push_back(dbe->book2D("Truth_PionMatch_dPtvslength"+tmp+amb,"Truth Pion Tau-Matchd dP_{t}^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,100,-5,5)); Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("P_{T}-P_{T}^{Truth} (GeV)",2); Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_PionMatch_dphivslength.at(ambiguity).push_back(dbe->book2D("Truth_PionMatch_dphivslength"+tmp+amb,"Truth Pion Tau-Matchd d#phi^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,64,-0.1,0.1)); Truth_PionMatch_dphivslength.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (rad)",2); Truth_PionMatch_dphivslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_PionMatch_dthetavslength.at(ambiguity).push_back(dbe->book2D("Truth_PionMatch_dthetavslength"+tmp+amb,"Truth Pion Tau-Matchd d#theta^{Lab} vs L( "+tmp+amb+")",100,0.0,5.0,64,-0.1,0.1)); Truth_PionMatch_dthetavslength.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{Truth} (rad)",2); Truth_PionMatch_dthetavslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_PionMatch_reldPtvslength.at(ambiguity).push_back(dbe->book2D("Truth_PionMatch_reldPtvslength"+tmp+amb,"Truth Pion Tau-Matchd Relative dP_{t}^{Lab}vs L( "+tmp+amb+")",100,0.0,5.0,100,-5,5)); Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_PionMatch_dPtvslength.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);


	  TruthVtxX.at(ambiguity).push_back(dbe->book1D("TruthPVtxX"+tmp+amb,"Vtx_{x,"+tmp+"}^{Prime,KF}-Vtx_{x}^{Prime,Truth} "+amb,100 ,-0.01,0.01)); axis="Vtx_{x,"+tmp+"}^{Prime,KF}-Vtx_{x}^{Prime,Truth} (cm)"; TruthVtxX.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthVtxY.at(ambiguity).push_back(dbe->book1D("TruthPVtxY"+tmp+amb,"Vtx_{y,"+tmp+"}^{Prime,KF}-Vtx_{y}^{Prime,Truth} "+amb,100 ,-0.01,0.01));axis="Vtx_{y,"+tmp+"}^{Prime,KF}-Vtx_{y}^{Prime,Truth} (cm)"; TruthVtxY.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthVtxZ.at(ambiguity).push_back(dbe->book1D("TruthPVtxZ"+tmp+amb,"Vtx_{z,"+tmp+"}^{Prime,KF}-Vtx_{z}^{Prime,Truth} "+amb,100 ,-0.01,0.01));axis="Vtx_{z,"+tmp+"}^{Prime,KF}-Vtx_{z}^{Prime,Truth} (cm)"; TruthVtxZ.at(ambiguity).at(idx)->setAxisTitle(axis.Data());


	  TruthSecVtxX.at(ambiguity).push_back(dbe->book1D("TruthSecVtxX"+tmp+amb,"Vtx_{x,"+tmp+"}^{Sec.,KF}-Vtx_{x}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{x,"+tmp+"}^{Sec.,KF}-Vtx_{x}^{Sec.,Truth} (cm)"; TruthSecVtxX.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxY.at(ambiguity).push_back(dbe->book1D("TruthSecVtxY"+tmp+amb,"Vtx_{y,"+tmp+"}^{Sec,KF}-Vtx_{y}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{y,"+tmp+"}^{Sec.,KF}-Vtx_{y}^{Sec.,Truth} (cm)"; TruthSecVtxY.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxZ.at(ambiguity).push_back(dbe->book1D("TruthSecVtxZ"+tmp+amb,"Vtx_{z,"+tmp+"}^{Sec,KF}-Vtx_{z}^{Sec.,Truth} "+amb,100 ,-0.1,0.1));axis="Vtx_{z,"+tmp+"}^{Sec.,KF}-Vtx_{z}^{Sec.,Truth} (cm)"; TruthSecVtxZ.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

          TruthSecVtxXp.at(ambiguity).push_back(dbe->book1D("TruthSecVtxXp"+tmp+amb,"Vtx_{x',"+tmp+"}^{Sec.,KF}-Vtx_{x'}^{Sec.,Truth} "+amb,100 ,-0.05,0.05));axis="Vtx_{x',"+tmp+"}^{Sec.,KF}-Vtx_{x'}^{Sec.,Truth} (cm)"; TruthSecVtxXp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxYp.at(ambiguity).push_back(dbe->book1D("TruthSecVtxYp"+tmp+amb,"Vtx_{y',"+tmp+"}^{Sec,KF}-Vtx_{y'}^{Sec.,Truth} "+amb,100 ,-0.05,0.05));axis="Vtx_{y',"+tmp+"}^{Sec.,KF}-Vtx_{y'}^{Sec.,Truth} (cm)"; TruthSecVtxYp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          TruthSecVtxZp.at(ambiguity).push_back(dbe->book1D("TruthSecVtxZp"+tmp+amb,"Vtx_{z',"+tmp+"}^{Sec,KF}-Vtx_{z'}^{Sec.,Truth} "+amb,100 ,-1.0,1.0));axis="Vtx_{z',"+tmp+"}^{Sec.,KF}-Vtx_{z'}^{Sec.,Truth} (cm)"; TruthSecVtxZp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

          PullVtxX.at(ambiguity).push_back(dbe->book1D("PullPVtxX"+tmp+amb,"Pull Vtx_{x,"+tmp+"}^{Prime,KF} "+amb,100 ,-5.0,5.0)); axis="Pull Vtx_{x,"+tmp+"}^{Prime,KF}"; PullVtxX.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullVtxY.at(ambiguity).push_back(dbe->book1D("PullPVtxY"+tmp+amb,"Pull Vtx_{y,"+tmp+"}^{Prime,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{y,"+tmp+"}^{Prime,KF}"; PullVtxY.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullVtxZ.at(ambiguity).push_back(dbe->book1D("PullPVtxZ"+tmp+amb,"Pull Vtx_{z,"+tmp+"}^{Prime,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{z,"+tmp+"}^{Prime,KF}"; PullVtxZ.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

          PullSecVtxX.at(ambiguity).push_back(dbe->book1D("PullSecVtxX"+tmp+amb,"Pull Vtx_{x,"+tmp+"}^{Sec.,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{x,"+tmp+"}^{Sec.,KF}"; PullSecVtxX.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullSecVtxY.at(ambiguity).push_back(dbe->book1D("PullSecVtxY"+tmp+amb,"Pull Vtx_{y,"+tmp+"}^{Sec,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{y,"+tmp+"}^{Sec.,KF}"; PullSecVtxY.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullSecVtxZ.at(ambiguity).push_back(dbe->book1D("PullSecVtxZ"+tmp+amb,"Pull Vtx_{z,"+tmp+"}^{Sec,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{z,"+tmp+"}^{Sec.,KF}"; PullSecVtxZ.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

          PullSecVtxXp.at(ambiguity).push_back(dbe->book1D("PullSecVtxXp"+tmp+amb,"Pull Vtx_{x',"+tmp+"}^{Sec.,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{x',"+tmp+"}^{Sec.,KF}"; PullSecVtxXp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullSecVtxYp.at(ambiguity).push_back(dbe->book1D("PullSecVtxYp"+tmp+amb,"Pull Vtx_{y',"+tmp+"}^{Sec,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{y',"+tmp+"}^{Sec.,KF}"; PullSecVtxYp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullSecVtxZp.at(ambiguity).push_back(dbe->book1D("PullSecVtxZp"+tmp+amb,"Pull Vtx_{z',"+tmp+"}^{Sec,KF} "+amb,100 ,-5.0,5.0));axis="Pull Vtx_{z',"+tmp+"}^{Sec.,KF}"; PullSecVtxZp.at(ambiguity).at(idx)->setAxisTitle(axis.Data());


          PullTauPx.at(ambiguity).push_back(dbe->book1D("PullTauPx"+tmp+amb,"Pull P_{x}^{#tau,"+tmp+"} "+amb,100 ,-5.0,5.0));axis="Pull P_{x}^{#tau}"; PullTauPx.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullTauPy.at(ambiguity).push_back(dbe->book1D("PullTauPy"+tmp+amb,"Pull P_{y}^{#tau,"+tmp+"} "+amb,100 ,-5.0,5.0));axis="Pull P_{y}^{#tau}"; PullTauPy.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
          PullTauPz.at(ambiguity).push_back(dbe->book1D("PullTauPz"+tmp+amb,"Pull P_{z}^{#tau,"+tmp+"} "+amb,100 ,-5.0,5.0));axis="Pull P_{z}^{#tau}"; PullTauPz.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

	  TruthPVtxSig.at(ambiguity).push_back(dbe->book1D("PVtxSig"+tmp+amb,"#sigma_{Prime Vtx,"+tmp+"}^{Truth} "+amb,100 ,0.0,10.0)); axis="#sigma_{Prime Vtx,"+tmp+"}^{Truth}"; TruthPVtxSig.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthSecVtxSig.at(ambiguity).push_back(dbe->book1D("SecVtxSig"+tmp+amb,"#sigma_{Sec. Vtx"+tmp+"}^{Truth} "+amb,100 ,0.0,10.0));axis="#sigma_{Sec. Vtx"+tmp+"}^{Truth} "; TruthSecVtxSig.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthTauFlightDir.at(ambiguity).push_back(dbe->book1D("TauFlightDir"+tmp+amb,"|d#psi_{#tau Dir.,#tau}^{KF}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| "+amb,100 ,0.0,0.1));axis="|d#psi_{#tau Dir.,#tau}^{KF}-d#psi_{#tau Dir.,#tau}^{Truth,"+tmp+"}| (rad)"; TruthTauFlightDir.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  TruthTauFlightDirCheck.at(ambiguity).push_back(dbe->book1D("TauFlightDirCheck"+tmp+amb,"|#psi_{#tau Dir.,Vtx}^{Truth}-#psi_{#tau Dir.,#tau}^{Truth}| "+amb,100 ,0.0,0.01));axis="|#psi_{#tau Dir.,Vtx}^{Truth}-#psi_{#tau Dir.,#tau}^{Truth}| (rad)"; TruthTauFlightDirCheck.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  
	  Truth_TauMatch_dPt.at(ambiguity).push_back(dbe->book1D("TaudResPt"+tmp+amb,"dP_{t,"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));   axis="dP_{t,"+tmp+"}^{#tau Match} (GeV)"; Truth_TauMatch_dPt.at(ambiguity).at(idx)->setAxisTitle(axis.Data());
	  Truth_TauMatch_dPz.at(ambiguity).push_back(dbe->book1D("TauResPz"+tmp+amb,"dP_{z,"+tmp+"}^{#tau Match} "+amb,100 ,-50.0,50.0));   axis="dP_{z,"+tmp+"}^{#tau Match} (GeV)"; Truth_TauMatch_dPz.at(ambiguity).at(idx)->setAxisTitle(axis.Data());

	  Truth_TauMatch_dPtvsL.at(ambiguity).push_back(dbe->book2D("TauResPtvsL"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",100,0.0,5.0,100,-50,50)); Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_TauMatch_dEvsL.at(ambiguity).push_back(dbe->book2D("TauEResvsL"+tmp+amb,"Truth Tau Matched E vs L "+tmp+amb+")",100,0.0,5.0,100,-50,50)); Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_TauMatch_dphivsL.at(ambiguity).push_back(dbe->book2D("TaudphivsL"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",100,0.0,5.0,100,-0.1,0.1)); Truth_TauMatch_dphivsL.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (rad)",2); Truth_TauMatch_dphivsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_TauMatch_dthetavsL.at(ambiguity).push_back(dbe->book2D("TauThetaResvsL"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",100,0.0,5.0,100,-0.1,0.1)); Truth_TauMatch_dthetavsL.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (rad)",2); Truth_TauMatch_dthetavsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
	  Truth_TauMatch_reldPtvsL.at(ambiguity).push_back(dbe->book2D("TauReldPtvsL"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs L( "+tmp+amb+")",100,0.0,5.0,100,-5,5)); Truth_TauMatch_reldPtvsL.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);

	  Truth_TauMatch_dGFAnglevsL.at(ambiguity).push_back(dbe->book2D("GFResAnglevsL"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs L( "+tmp+amb+")",10,0.0,5.0,64,-0.1,0.1)); Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsL.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);

	  Truth_TauMatch_dGFAngle.at(ambiguity).push_back(dbe->book1D("GFResAngle"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} "+tmp+amb+")",100,-0.025,0.025)); Truth_TauMatch_dGFAngle.at(ambiguity).at(idx)->setAxisTitle("d#theta_{GF,#tau} (rad)"); 

	  Truth_TauMatch_TrueMaxGFAngle.at(ambiguity).push_back(dbe->book1D("GFMaxAngle"+tmp+amb,"Truth Tau-Matched #theta_{GF}^{Max,Truth} "+tmp+amb+")",100,0.0,0.1)); Truth_TauMatch_TrueMaxGFAngle.at(ambiguity).at(idx)->setAxisTitle("#theta_{GF}^{Max,Truth} (rad)");
	  Truth_TauMatch_TrueGFAngle.at(ambiguity).push_back(dbe->book1D("GFAngle"+tmp+amb,"Truth Tau-Matched #theta_{GF}^{Truth} "+tmp+amb+")",100,0.0,0.1)); Truth_TauMatch_TrueGFAngle.at(ambiguity).at(idx)->setAxisTitle("#theta_{GF}^{Truth} (rad)");

	  Truth_TauMatch_TrueGFAngle_TauRF.at(ambiguity).push_back(dbe->book1D("GFAngle_TauRF"+tmp+amb,"Truth Tau-Matched cos(#theta_{GF}^{Truth}) "+tmp+amb+")",100,-1.0,1.0)); Truth_TauMatch_TrueGFAngle_TauRF.at(ambiguity).at(idx)->setAxisTitle("cos(#theta_{GF,#tau RF}^{Truth})");

	  Truth_TauMatch_GFAngleRecovsTruth_TauRF.at(ambiguity).push_back(dbe->book2D("GFAngleRecovsTruth_TauRF"+tmp+amb,"Truth Tau-Matched cos(#theta_{GF,#tau RF}^{Truth}) vs cos(#theta_{GF,#tau RF})",50,-1.0,1.0,50,-1.0,1.0)); Truth_TauMatch_GFAngleRecovsTruth_TauRF.at(ambiguity).at(idx)->setAxisTitle("cos(#theta_{GF,#tau RF}^{Truth})",2); Truth_TauMatch_GFAngleRecovsTruth_TauRF.at(ambiguity).at(idx)->setAxisTitle("cos(#theta_{GF,#tau RF})",1);


	  Truth_TauMatch_TrueGFAngleoverMaxGFAngle.at(ambiguity).push_back(dbe->book1D("GFAngleoverMaxGFAngle"+tmp+amb,"Truth Tau-Matched #theta_{GF}^{Truth}/#theta_{GF}^{Max,Truth} "+tmp+amb+")",100,0.0,1.0)); Truth_TauMatch_TrueGFAngleoverMaxGFAngle.at(ambiguity).at(idx)->setAxisTitle("#theta_{GF}^{Truth}/#theta_{GF}^{Max,Truth}");

	  Truth_TauMatch_dGFAnglevsTrackchi2.at(ambiguity).push_back(dbe->book2D("GFResAnglevsTrackchi2"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs prob( "+tmp+amb+")",20,0.0,1.0,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsTrackchi2.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsTrackchi2.at(ambiguity).at(idx)->setAxisTitle("Track Probability",1);
	  Truth_TauMatch_dGFAnglevsTrackQuality.at(ambiguity).push_back(dbe->book2D("GFResAnglevsTrackQuality"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs quality( "+tmp+amb+")",5,-0.5,4.5,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsTrackQuality.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsTrackQuality.at(ambiguity).at(idx)->setAxisTitle("Track Quaility",1);
	  Truth_TauMatch_dGFAnglevsMaxTrackPt.at(ambiguity).push_back(dbe->book2D("GFResAnglevsMaxTrackPt"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs P_{t}^{max}( "+tmp+amb+")",100,0.0,25.0,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsMaxTrackPt.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsMaxTrackPt.at(ambiguity).at(idx)->setAxisTitle("P_{t}^{max} (GeV)",1);
	  Truth_TauMatch_dGFAnglevsMinTrackPt.at(ambiguity).push_back(dbe->book2D("GFResAnglevsMinTrackPt"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs P_{t}^{min}( "+tmp+amb+")",40,0.0,10.0,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsMinTrackPt.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsMinTrackPt.at(ambiguity).at(idx)->setAxisTitle("P_{t}^{min} (GeV)",1);
	  Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt.at(ambiguity).push_back(dbe->book2D("GFResAnglevsMinoverMaxTrackPt"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs P_{t}^{min}/P_{t}^{max}( "+tmp+amb+")",20,0.0,1.0,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsMinoverMaxTrackPt.at(ambiguity).at(idx)->setAxisTitle("P_{t}^{min}/P_{t}^{max}",1);
	  Truth_TauMatch_dGFAnglevsHitFrac.at(ambiguity).push_back(dbe->book2D("GFResAnglevsHitFrac"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs Hit fraction( "+tmp+amb+")",22,0.0,1.1,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsHitFrac.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsHitFrac.at(ambiguity).at(idx)->setAxisTitle("HitFrac",1);
	  Truth_TauMatch_dGFAnglevsVertexchi2.at(ambiguity).push_back(dbe->book2D("GFResAnglevsVertexchi2"+tmp+amb,"Truth Tau-Matched d#theta_{GF}^{Lab} vs prob( "+tmp+amb+")",20,0.0,1.0,64,-0.05,0.05)); Truth_TauMatch_dGFAnglevsVertexchi2.at(ambiguity).at(idx)->setAxisTitle("d#theta_{t,#tau} (rad)",2); Truth_TauMatch_dGFAnglevsVertexchi2.at(ambiguity).at(idx)->setAxisTitle("Vertex Probability",1);


	  Truth_TauMatch_dPtvsTheta.at(ambiguity).push_back(dbe->book2D("TauResPtvsTheta"+tmp+amb,"Truth Tau-Matched P_{t} vs #theta( "+tmp+amb+")",32,0.0,TMath::Pi(),100,-50,50)); Truth_TauMatch_dPtvsTheta.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta (rad)",1);
          Truth_TauMatch_dEvsTheta.at(ambiguity).push_back(dbe->book2D("TauEResvsTheta"+tmp+amb,"Truth Tau Matched E vs #theta "+tmp+amb+")",32,0.0,TMath::Pi(),100,-50,50)); Truth_TauMatch_dEvsTheta.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta (rad)",1);
          Truth_TauMatch_dphivsTheta.at(ambiguity).push_back(dbe->book2D("TaudphivsTheta"+tmp+amb,"Truth Tau-Matched P_{t} vs #theta( "+tmp+amb+")",32,0.0,TMath::Pi(),100,-0.1,0.1)); Truth_TauMatch_dphivsTheta.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (rad)",2); Truth_TauMatch_dphivsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta (rad)",1);
          Truth_TauMatch_dthetavsTheta.at(ambiguity).push_back(dbe->book2D("TauThetaResvsTheta"+tmp+amb,"Truth Tau-Matched P_{t} vs #theta( "+tmp+amb+")",32,0.0,TMath::Pi(),100,-0.1,0.1)); Truth_TauMatch_dthetavsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (rad)",2); Truth_TauMatch_dthetavsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta (rad)",1);
          Truth_TauMatch_reldPtvsTheta.at(ambiguity).push_back(dbe->book2D("TauReldPtvsTheta"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs #theta("+tmp+amb+")",32,0.0,TMath::Pi(),100,-5,5)); Truth_TauMatch_reldPtvsTheta.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsTheta.at(ambiguity).at(idx)->setAxisTitle("#theta (rad)",1);

          Truth_TauMatch_dPtvsEta.at(ambiguity).push_back(dbe->book2D("TauResPtvsEta"+tmp+amb,"Truth Tau-Matched P_{t} vs #eta( "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dPtvsEta.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsEta.at(ambiguity).at(idx)->setAxisTitle("#eta (rad)",1);
          Truth_TauMatch_dEvsEta.at(ambiguity).push_back(dbe->book2D("TauEResvsEta"+tmp+amb,"Truth Tau Matched E vs #eta "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dEvsEta.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsEta.at(ambiguity).at(idx)->setAxisTitle("#eta (rad)",1);
          Truth_TauMatch_dphivsEta.at(ambiguity).push_back(dbe->book2D("TaudphivsEta"+tmp+amb,"Truth Tau-Matched P_{t} vs #eta( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dphivsEta.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (rad)",2); Truth_TauMatch_dphivsEta.at(ambiguity).at(idx)->setAxisTitle("#eta (rad)",1);
          Truth_TauMatch_dthetavsEta.at(ambiguity).push_back(dbe->book2D("TauThetaResvsEta"+tmp+amb,"Truth Tau-Matched P_{t} vs #eta( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dthetavsEta.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (rad)",2); Truth_TauMatch_dthetavsEta.at(ambiguity).at(idx)->setAxisTitle("#eta (rad)",1);
          Truth_TauMatch_reldPtvsEta.at(ambiguity).push_back(dbe->book2D("TauReldPtvsEta"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs #eta("+tmp+amb+")",32,-2,2,100,-1,1)); Truth_TauMatch_reldPtvsEta.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsEta.at(ambiguity).at(idx)->setAxisTitle("#eta (rad)",1);

	  ////////////////////// Pt
          Truth_TauMatch_dPtvsPt.at(ambiguity).push_back(dbe->book2D("TauResPtvsPt"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dPtvsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dEvsPt.at(ambiguity).push_back(dbe->book2D("TauEResvsPt"+tmp+amb,"Truth Tau Matched E vs P_{t} "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dEvsPt.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dphivsPt.at(ambiguity).push_back(dbe->book2D("TaudphivsPt"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dphivsPt.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (GeV)",2); Truth_TauMatch_dphivsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dthetavsPt.at(ambiguity).push_back(dbe->book2D("TauThetaResvsPt"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dthetavsPt.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (GeV)",2); Truth_TauMatch_dthetavsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_reldPtvsPt.at(ambiguity).push_back(dbe->book2D("TauReldPtvsPt"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs P_{t}("+tmp+amb+")",32,-2,2,100,-1,1)); Truth_TauMatch_reldPtvsPt.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsPt.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);

	  ///////////////////////// neta<1.0
          Truth_TauMatch_dPtvsPtandEtaCut.at(ambiguity).push_back(dbe->book2D("TauResPtvsPtandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dPtvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dEvsPtandEtaCut.at(ambiguity).push_back(dbe->book2D("TauEResvsPtandEtaCut"+tmp+amb,"Truth Tau Matched E vs P_{t} "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dEvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dphivsPtandEtaCut.at(ambiguity).push_back(dbe->book2D("TaudphivsPtandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dphivsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (GeV)",2); Truth_TauMatch_dphivsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_dthetavsPtandEtaCut.at(ambiguity).push_back(dbe->book2D("TauThetaResvsPtandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs P_{t}( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dthetavsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (GeV)",2); Truth_TauMatch_dthetavsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);
          Truth_TauMatch_reldPtvsPtandEtaCut.at(ambiguity).push_back(dbe->book2D("TauReldPtvsPtandEtaCut"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs P_{t}("+tmp+amb+")",32,-2,2,100,-1,1)); Truth_TauMatch_reldPtvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsPtandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t} (GeV)",1);

          ///////////////////////// neta<1.0
          Truth_TauMatch_dPtvsLandEtaCut.at(ambiguity).push_back(dbe->book2D("TauResPtvsLandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dPtvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
          Truth_TauMatch_dEvsLandEtaCut.at(ambiguity).push_back(dbe->book2D("TauEResvsLandEtaCut"+tmp+amb,"Truth Tau Matched E vs L "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dEvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
          Truth_TauMatch_dphivsLandEtaCut.at(ambiguity).push_back(dbe->book2D("TaudphivsLandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dphivsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (GeV)",2); Truth_TauMatch_dphivsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
          Truth_TauMatch_dthetavsLandEtaCut.at(ambiguity).push_back(dbe->book2D("TauThetaResvsLandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t} vs L( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dthetavsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (GeV)",2); Truth_TauMatch_dthetavsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);
          Truth_TauMatch_reldPtvsLandEtaCut.at(ambiguity).push_back(dbe->book2D("TauReldPtvsLandEtaCut"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab}vs L("+tmp+amb+")",32,-2,2,100,-1,1)); Truth_TauMatch_reldPtvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsLandEtaCut.at(ambiguity).at(idx)->setAxisTitle("L (cm)",1);


	  ///////////////////////// neta<1.0
          Truth_TauMatch_dPtvsPhiandEtaCut.at(ambiguity).push_back(dbe->book2D("TauResPtvsPhiandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t}  vs #phi( "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dPtvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{t,#tau} (GeV)",2); Truth_TauMatch_dPtvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi (rad)",1);
          Truth_TauMatch_dEvsPhiandEtaCut.at(ambiguity).push_back(dbe->book2D("TauEResvsPhiandEtaCut"+tmp+amb,"Truth Tau Matched E  vs #phi "+tmp+amb+")",32,-2,2,100,-50,50)); Truth_TauMatch_dEvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("E_{#tau} (GeV)",2); Truth_TauMatch_dEvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi (rad)",1);
          Truth_TauMatch_dphivsPhiandEtaCut.at(ambiguity).push_back(dbe->book2D("TaudphivsPhiandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t}  vs #phi( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dphivsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi-#phi^{Truth} (GeV)",2); Truth_TauMatch_dphivsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi (rad)",1);
          Truth_TauMatch_dthetavsPhiandEtaCut.at(ambiguity).push_back(dbe->book2D("TauThetaResvsPhiandEtaCut"+tmp+amb,"Truth Tau-Matched P_{t}  vs #phi( "+tmp+amb+")",32,-2,2,100,-0.1,0.1)); Truth_TauMatch_dthetavsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#theta-#theta^{truth} (GeV)",2); Truth_TauMatch_dthetavsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi (rad)",1);
          Truth_TauMatch_reldPtvsPhiandEtaCut.at(ambiguity).push_back(dbe->book2D("TauReldPtvsPhiandEtaCut"+tmp+amb,"Truth Tau Tau-Matchd Relative dP_{t}^{Lab} vs #phi("+tmp+amb+")",32,-2,2,100,-1,1)); Truth_TauMatch_reldPtvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("P_{T}/P_{T}^{Truth}-1 (GeV)",2); Truth_TauMatch_reldPtvsPhiandEtaCut.at(ambiguity).at(idx)->setAxisTitle("#phi (rad)",1);


	  Truth_PionDr.at(ambiguity).push_back(dbe->book1D("Truth_PionDr"+amb,"Truth Pion Dr",100, 0.0,0.5));  Truth_PionDr.at(ambiguity).at(idx)->setAxisTitle("Pion Dr_{max}^{Truth}");

	  Truth_PionDrAll.at(ambiguity).push_back(dbe->book1D("Truth_PionDrAll"+amb,"Truth Pion Dr",100, 0.0,0.5));  Truth_PionDrAll.at(ambiguity).at(idx)->setAxisTitle("Pion Dr_{max}^{Truth,All}");

	  Truth_a1MassHPS.at(ambiguity).push_back(dbe->book1D("Truth_a1MassHPS"+amb,"Truth_a1MassHPS",100, 0.5,0.5));  Truth_a1MassHPS.at(ambiguity).at(idx)->setAxisTitle("M_{a1}^{HPS}");

          Truth_MatchPionDr.at(ambiguity).push_back(dbe->book1D("Truth_MatchPionDr"+amb,"Truth Pion Dr",100, 0.0,0.5));        Truth_MatchPionDr.at(ambiguity).at(idx)->setAxisTitle("Pion Dr_{max}^{Truth}");
          Truth_MatchPionDrHPS.at(ambiguity).push_back(dbe->book1D("Truth_MatchPionDrHPS"+amb,"Truth Pion Dr",100, 0.0,0.5));  Truth_MatchPionDrHPS.at(ambiguity).at(idx)->setAxisTitle("Pion Dr_{max}^{Truth,All}");



	} 
      }
    }
  }
}

void KinematicTauAnalyzer::endJob(){
  for(unsigned int i=0;i<MultiProngTauSolver::NAmbiguity;i++){
    float ratio = 0.0;
    if(cnt_!=0) ratio=(float)cntFound_.at(i)/cnt_;
    std::cout << "KinematicTauAnalyzer Ambiguity " << i <<"-->  Efficiency: "<< cntFound_.at(i)<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%" << std::endl;
  }
}
bool KinematicTauAnalyzer::doJAKID(unsigned int i){
  for(unsigned int j=0;j<JAKID_.size();j++){
    if(((unsigned int)JAKID_.at(j))==i)return true;
  }
  return false;
}

bool KinematicTauAnalyzer::A1Matched(std::vector<const reco::GenParticle* > &DecayProd,TLorentzVector a1){
  TLorentzVector mc_vis(0,0,0,0);
  for(unsigned int j=0;j<DecayProd.size();j++){
    if(fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::tau_plus)){
      TLorentzVector LV(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
      mc_vis+=LV;
    }
    if(fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::nu_e) || fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::nu_mu) || fabs(DecayProd.at(j)->pdgId())==fabs(PdtPdgMini::nu_tau) ){
      TLorentzVector LV(DecayProd.at(j)->p4().Px(),DecayProd.at(j)->p4().Py(),DecayProd.at(j)->p4().Pz(),DecayProd.at(j)->p4().E());
      mc_vis-=LV;
    }
  }
  return a1.DeltaR(mc_vis)<TauMatchingDR_;
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauAnalyzer);

