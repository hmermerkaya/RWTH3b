#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"

int ThreeProngTauCreator::create(unsigned int &ambiguity,SelectedKinematicDecay &KFTau){
  status=0;
  std::cout << "Fit start" << std::endl;
  std::vector<TrackParticle> pions;
  std::vector<LorentzVectorParticle> neutrino;
  std::vector<LorentzVectorParticle> daughters;
  // Helix fit for a1
  ConfigurePions(KFTau,pions);
  if(pions.size()!=3)return status;
  SecondaryVertexHelper SVH(transientTrackBuilder_,KFTau);
  TransientVertex secVtx=SVH.InitialSecondaryVertex();
  std::cout << "Fit start before vertex fit" << std::endl;
  if(!FitA1(pions,/*secVtx*/PV_)){return status;}
  std::cout << "Fit start after vertex fit" << std::endl;
  A1=mother();
  // Tau fit
  daughters.push_back(A1);
  std::cout << "Fit start before Tau fit" << std::endl;
  if(!FitTau(daughters,PV_,ambiguity)){return status;}
  std::cout << "Fit start After vertex fit" << std::endl;
  Tau=mother();
  // Fit sequence complete
  status=1;
  std::cout << " Fit Complete" << std::endl;
  return status;
}

void ThreeProngTauCreator::ConfigurePions(SelectedKinematicDecay &KFTau, std::vector<TrackParticle> &pions){
  PV_=KFTau.InitialPrimaryVertexReFit();
  GlobalPoint pv(PV_.position().x(),PV_.position().y(),PV_.position().z());
  std::vector<reco::TrackRef> selectedTracks=KFTau.InitialTrackTriplet();
  for(unsigned int i = 0; i!=selectedTracks.size();i++){
    pions.push_back(ParticleBuilder::CreateTrackParticle(selectedTracks.at(i),transientTrackBuilder_,pv));
  }
}

bool ThreeProngTauCreator::FitTau(std::vector<LorentzVectorParticle>  &unfitDaughters,const reco::Vertex &primaryVertex,
				  unsigned int &ambiguity){
  if(unfitDaughters.size()!=1){
    edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
    return false;
  }
  // Setup Constraint
  TVector3 pv(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);
  math::Error<LorentzVectorParticle::NVertex>::type pvCov;
  primaryVertex.fill(pvCov);
  for(int i = 0; i<LorentzVectorParticle::NVertex; i++){
    for(int j = 0; j<LorentzVectorParticle::NVertex; j++){
      pvcov(i,j)=pvCov(i,j);
      pvcov(j,i)=pvCov(i,j);
    }
  }
  TauA1NuConstrainedFitter TauA1NU(ambiguity,unfitDaughters.at(0),pv,pvcov,PMH.Get_tauMass());
  TauA1NU.SetMaxDelta(0.01);
  TauA1NU.SetNIterMax(1000);
  std::cout << "before fit" << std::endl;
  TauA1NU.Fit();
  std::cout << "after fit" << std::endl;
  if (TauA1NU.isConverged()) {
    //////////////////////////////////////////////////////////////////////////////////
    // debug statments 
    /*
    std::vector<LorentzVectorParticle> fitDaughters=TauA1NU.GetReFitDaughters();
    LorentzVectorParticle Mother=TauA1NU.GetMother();

    for(reco::GenParticleCollection::const_iterator itr = GenPart_->begin(); itr!= GenPart_->end(); ++itr){
      if(fabs(itr->pdgId())==15){
	const reco::GenParticle mytau=(*itr);
	for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){
	  const reco::Candidate *dau=(itr)->daughter(i);
	  if(fabs(dau->pdgId())==20213){
	    std::cout << "\n ============================" << std::endl;
	    std::cout << "Tau " << itr->p4().E() << " " << (itr)->vx() << " " << (itr)->vy() << " " << (itr)->vz()
                      << " "    << itr->p4().Px() << " " << itr->p4().Py() << " " << itr->p4().Pz() << " " << itr->p4().M() 
		      << " theta " << itr->p4().Theta() << " phi " << itr->p4().Phi()  <<  std::endl;
	    std::cout << "    ";
	    for(int i=-1; i<LorentzVectorParticle::NLorentzandVertexPar;i++){std::cout <<Mother.Parameter(i) << " ";}
	    std::cout << "\n ============================" << std::endl;
	    std::cout << "A1  " << dau->p4().E() << " " <<  dau->vx()  << " " << dau->vy() << " " << dau->vz() 
		      << " "    << dau->p4().Px() << " " << dau->p4().Py() << " " << dau->p4().Pz() << " " << dau->p4().M() <<  std::endl;
	    std::cout << "    ";
	    for(int k=-1; k<LorentzVectorParticle::NLorentzandVertexPar;k++){std::cout << fitDaughters.at(0).Parameter(k) << " ";}
	    std::cout << "\n ============================" << std::endl;
	  }
	  if(fabs(dau->pdgId())==16){
	    std::cout << "\n ============================" << std::endl;
	    std::cout << "Nu  " << dau->p4().E() << " " <<  dau->vx()  << " " << dau->vy() << " " << dau->vz()
                      << " "    << dau->p4().Px() << " " << dau->p4().Py() << " " << dau->p4().Pz() << " " << dau->p4().M() <<  std::endl;
	    std::cout << "    ";
            for(int k=-1; k<LorentzVectorParticle::NLorentzandVertexPar;k++){std::cout << fitDaughters.at(1).Parameter(k) << " ";}
	    std::cout << "\n ============================" << std::endl;
	  }
	}
      }
    }
    
    ///////////////////////////////////////////////////
    for(unsigned int i=0;i<fitDaughters.size();i++){
      for(unsigned int j=0;j<unfitDaughters.size();j++){
        fitDaughters.erase (fitDaughters.begin()+i);
        fitDaughters.push_back(unfitDaughters.at(j));
      }
    }
    */
    /////////////////////
    StoreResults(TauA1NU.ChiSquare(),TauA1NU.NDF(),TauA1NU.CSum(),TauA1NU.NIter(),TauA1NU.NConstraints(),TauA1NU.GetReFitDaughters()/*fitDaughters*/,TauA1NU.GetMother());
    LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Valid tree.";
    return true;
  } 
  LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: Warning! Tree is not converged. Skip tauCand.";//DEBUG
  return false;
}

bool ThreeProngTauCreator::FitA1(std::vector<TrackParticle> &pions,const reco::Vertex & primaryVertex){
  TVector3 pv(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  Chi2VertexFitter chi2v(pions,pv);
  chi2v.Fit();
  double c(0); for(unsigned int i=0;i<pions.size();i++){c+=pions.at(i).Charge();}
  int pdgid=fabs(PdtPdgMini::a_1_plus)*c;

  ////////////////////////////////////////////////////////////
  // debug
  /*
  LorentzVectorParticle LV=chi2v.GetMother(pdgid);
  TLorentzVector A1;
  double vx,vy,vz;
  for(reco::GenParticleCollection::const_iterator itr = GenPart_->begin(); itr!= GenPart_->end(); ++itr){
    if(fabs(itr->pdgId())==15){
      const reco::GenParticle mytau=(*itr);
      for (unsigned int i=0; i<(itr)->numberOfDaughters();i++){
	const reco::Candidate *dau=(itr)->daughter(i);
	if(fabs(dau->pdgId())==20213){
	  TLorentzVector a1(dau->p4().Px(),dau->p4().Py(),dau->p4().Pz(),dau->p4().E());
	  A1=a1;
	  //if(a1.DeltaR(LV)<0.4){
	    // Print info from matching tau
	  std::cout << "A1  Px " << dau->p4().Px() << " Py " << dau->p4().Py() << " Pz " << dau->p4().Pz() << " E " << dau->p4().E() << " M " << dau->p4().M() <<  std::endl;
	  std::cout << "A1  vx "       << dau->vx()  << " vy " << dau->vy() << " vz " << dau->vz() << std::endl;
	  std::cout << "Tau Px " << itr->p4().Px() << " Py " << itr->p4().Py() << " Pz " << itr->p4().Pz() << " E " << itr->p4().E() << " M " << itr->p4().M() <<  std::endl;
	  std::cout << "Tau vx "      << (itr)->vx() << " vy " << (itr)->vy() << " vz " << (itr)->vz() << std::endl;
	  vx=(dau)->vx();
	  vy=(dau)->vy();
	  vz=(dau)->vz();

	  reco::Vertex::Error ve=PV_.covariance() ;
	  reco::Vertex::Point vp((itr)->vx(),(itr)->vy(),(itr)->vz());
	  PV_=reco::Vertex(vp,ve);
	  //}
	}
      }
    }
  }
  std::cout << "Chi2Vertex Fit vx " << chi2v.GetVertex().X() << " vy " << chi2v.GetVertex().Y() << " " << chi2v.GetVertex().Z() << " " << std::endl; 
  //std::cout << "Chi2Vertex a1 Px  " << LV.Px() << " py " << LV.Py() << " pz " << LV.Pz() << std::endl; 
  TMatrixT<double> par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> cov(LorentzVectorParticle::NLorentzandVertexPar);

  for(unsigned int i=0; i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    par(i,0)=LV.Parameter(i);
    for(unsigned int j=0; j<LorentzVectorParticle::NLorentzandVertexPar;j++){cov(i,j)=LV.Covariance(i,j);}
  }
  par(LorentzVectorParticle::vx,0)=vx;
  par(LorentzVectorParticle::vy,0)=vy;
  par(LorentzVectorParticle::vz,0)=vz;
  par(LorentzVectorParticle::px,0)=A1.Px();
  par(LorentzVectorParticle::py,0)=A1.Py();
  par(LorentzVectorParticle::pz,0)=A1.Pz();
  par(LorentzVectorParticle::m,0)=A1.M();
  LorentzVectorParticle newA1(par,cov, LV.PDGID(),LV.Charge(), LV.BField());
  */
  ////////////////////////////////////////////////////////////
  StoreResults(chi2v.ChiSquare(),chi2v.NDF(),0,0,0,chi2v.GetReFitLorentzVectorParticles(),chi2v.GetMother(pdgid));  
  return  true;
}

