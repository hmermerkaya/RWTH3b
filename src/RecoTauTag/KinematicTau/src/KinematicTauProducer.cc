#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"


KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
  fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
  primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag")),
  trkCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauDaughterTracks" ) ),
  cntSVFound_(MultiProngTauSolver::NAmbiguity),
  cntSVQC_(MultiProngTauSolver::NAmbiguity),
  cntLCFit_(MultiProngTauSolver::NAmbiguity),
  cntLCFitQC_(MultiProngTauSolver::NAmbiguity),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  minTau_(iConfig.getUntrackedParameter<unsigned int>("minTau", 1)),
  etacut_(iConfig.getUntrackedParameter<double>("etacut",2.1)),
  sigcut_(iConfig.getUntrackedParameter<double>("sigcut",3.0)),
  do_BDTTrain_(iConfig.getUntrackedParameter("do_BDTTrain",(bool)(false))),
  do_BDTComp_(iConfig.getUntrackedParameter("do_BDTComp",(bool)(true))),
  BDTweightFileMinus_(iConfig.getUntrackedParameter<std::string>("BDTweightFileMinus")),
  BDTweightFilePlus_(iConfig.getUntrackedParameter<std::string>("BDTweightFilePlus")),
  BDTweightFileZero_(iConfig.getUntrackedParameter<std::string>("BDTweightFileZero"))
{
  produces<reco::RecoChargedCandidateCollection>("KinematicFitTauDaughters");
  produces<SelectedKinematicDecayCollection>("KinematicFitTau");
}

KinematicTauProducer::~KinematicTauProducer(){
}

// ------------ method called on each new Event  ------------
void KinematicTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool filterValue = false;
  cnt_++;
  iEvent_ = &iEvent;
  std::auto_ptr<SelectedKinematicDecayCollection> KinematicFitTauDecays = std::auto_ptr<SelectedKinematicDecayCollection >(new SelectedKinematicDecayCollection);
  SelectedKinematicDecayCollection &KinematicFitTauDecays_=*KinematicFitTauDecays;
  filterValue = select(KinematicFitTauDecays_,iSetup);
  if(filterValue) edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Passed";
  else edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Failed";
  iEvent_->put(KinematicFitTauDecays,"KinematicFitTau");
}

void KinematicTauProducer::beginJob(){
  cnt_ = 0;
  cntFound_=0;
  if(do_BDTTrain_){
   output = new TFile("ForTrain.root","RECREATE");
   output_tree = new TTree("t","t");
   output_tree->Branch("BDT_vtxSignPVRotSV",&BDT_vtxSignPVRotSV);
   output_tree->Branch("BDT_vtxSignPVRotPVRed",&BDT_vtxSignPVRotPVRed);
   output_tree->Branch("BDT_a1Mass",&BDT_a1Mass);
   output_tree->Branch("BDT_energyTFraction",&BDT_energyTFraction);
   output_tree->Branch("BDT_chiSquared",&BDT_chiSquared);
  }
  if(do_BDTComp_){
     reader  = new TMVA::Reader();
     reader->AddVariable( "fracMins", &fracMins);
     reader->AddVariable( "a1MassMins", &a1MassMins);
     reader->AddVariable( "ProbMins3", &ProbMins3);
     reader->AddVariable( "iterMins", &iterMins);
     reader->AddVariable( "PVSVMins", &PVSVMins);
     reader->BookMVA( "BDT_method",BDTweightFileZero_);
  }
  //     reader->BookMVA( "BDT_method", BDTweightFileMinus_ );
}

void KinematicTauProducer::endJob(){
  edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] Found 3 tracks Candidate Selection efficiency: "
				  << cntFound_ <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)cntFound_/(float)cnt_)*100.0
				  << "+/-" <<  sqrt((float)cntFound_*(1-(float)cntFound_/(float)cnt_))/(float)cnt_ <<"%";
  for(unsigned int i=0;i<MultiProngTauSolver::NAmbiguity;i++){
    std::cout << "===========================================================================================================" << std::endl; 
    edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] Found Secondary Vertex from 3 tracks Selection efficiency: "
				    << cntSVFound_.at(i) <<"/"<<cnt_<<" = "<<std::setprecision(4)<< (float)cntSVFound_.at(i)/(float)cnt_*100.0
				    << "+/-" <<  sqrt((float)cntSVFound_.at(i)*(1-(float)cntSVFound_.at(i)/(float)cnt_))/(float)cnt_ <<"%";
    edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] Found Secondary Vertex QC Selection efficiency: "
                                    << cntSVQC_.at(i) <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)cntSVQC_.at(i)/(float)cnt_)*100.0
				    << "+/-" <<  sqrt((float)cntSVQC_.at(i)*(1-(float)cntSVQC_.at(i)/(float)cnt_))/(float)cnt_ <<"%";
    edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] Found LC Fit succeeded Selection efficiency: "
                                    << cntLCFit_.at(i) <<"/"<<cnt_<<" = "<<std::setprecision(4)<< (float)cntLCFit_.at(i)/(float)cnt_*100.0
				    << "+/-" <<  sqrt((float)cntLCFit_.at(i)*(1-(float)cntLCFit_.at(i)/(float)cnt_))/(float)cnt_ <<"%";
    edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] Found LC Fit succeeded Selection efficiency: "
                                    << cntLCFitQC_.at(i) <<"/"<<cnt_<<" = "<<std::setprecision(4)<< (float)cntLCFitQC_.at(i)/(float)cnt_*100.0
				    << "+/-" <<  sqrt((float)cntLCFitQC_.at(i)*(1-(float)cntLCFitQC_.at(i)/(float)cnt_))/(float)cnt_ <<"%";

  }
  std::cout << "===========================================================================================================" << std::endl;
  if(do_BDTTrain_){
    output->Write();
    output->Close();
  }

}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,const edm::EventSetup& iSetup){
  //std::cout << "KinematicTauProducer::select" << std::endl;
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent_->getByLabel(gensrc_, genParticles);

  // Setup vertex for looking at Taus
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtxTag_, primaryVertexCollection);
  if (!primaryVertexCollection.isValid()) return false;
  reco::VertexCollection primaryVertices = *primaryVertexCollection;
  if(primaryVertices.size()>0){ // if event has vertex get tau candidated
    //std::cout << "found vertex" << std::endl;
    edm::Handle<std::vector<std::vector<SelectedKinematicDecay> > > KinematicTauCandidate;
    iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidate);
    //std::cout << "N candidates List " << KinematicTauCandidate->size() << std::endl;
    for(unsigned int i=0; i<KinematicTauCandidate->size();i++){
      std::vector<SelectedKinematicDecay>  KFTauCandidates;
      //std::cout << "Kinematic Fit " << i << std::endl;
      //std::cout << "N candidates" << KinematicTauCandidate->at(i).size() << std::endl;
      for(unsigned int j=0; j<KinematicTauCandidate->at(i).size();j++){
	//std::cout << "Kinematic Fit " << i << " " << j << std::endl;
	SelectedKinematicDecay KTau=KinematicTauCandidate->at(i).at(j);
        SecondaryVertexHelper SVH(transTrackBuilder_,KTau);
        if(SVH.hasSecondaryVertex()){
	  //////////////////////////////////////////////////////
	  //Rebuild primary vertex
	  reco::TrackCollection NonTauTracksLists_;
	  GetNonTauTracksFromVertex(KinematicTauCandidate->at(i).at(j),trkCollectionTag_,NonTauTracksLists_);
	  TransientVertex tmpVtx_;
	  std::vector<reco::TransientTrack> trks_;
	  for (reco::TrackCollection::iterator iter=NonTauTracksLists_.begin(); iter!=NonTauTracksLists_.end(); ++iter){
	    trks_.push_back(transTrackBuilder_->build(*iter));
	  }
	  if (!SecondaryVertexHelper::checkSecVtx(trks_,tmpVtx_,true))continue;
	  reco::Vertex primaryVertexReFit=tmpVtx_;
	  //////////////////////////////////////////////////////
	  reco::Vertex primaryVertexReFitAndRotated=primaryVertexReFit;
          TVector3 tauFlghtDirNoCorr;
          TVector3 tauFlghtDir;
          TLorentzVector a1_p4=SVH.Initial_a1_p4();
          double initThetaGJ,ThetaMax;
          TransientVertex SecondaryVertex=SVH.SecondaryVertex();
	  std::vector<reco::TransientTrack> RefittedTracks=SVH.RefittedTracks();
          double s = VertexRotationAndSignificance(SecondaryVertex,RefittedTracks,tauFlghtDirNoCorr,primaryVertexReFitAndRotated,a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	  if(/*(s<sigcut_  || fabs(initThetaGJ)<fabs(ThetaMax)) &&*/ s>=0 ){
	    if(fabs(tauFlghtDirNoCorr.Eta())<etacut_ || fabs(tauFlghtDir.Eta())<etacut_ ){
	      KTau.SetInitialVertexProperties(primaryVertexReFit,primaryVertexReFitAndRotated,SVH.RefittedTracks(),SVH.SecondaryVertex());
	      KTau.SetInitialKinematics(tauFlghtDirNoCorr,SVH.Initial_pions(),a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	      //std::cout << "before Fit "<< std::endl;
	      cntFound_++;
	      bool FitOK=FitKinematicTauCandidate(KTau,transTrackBuilder_,genParticles);
	      //std::cout << "after Fit "<< std::endl;
	      if(FitOK){
		//std::cout << "good Fit "<< std::endl;
		KFTauCandidates.push_back(KTau);
	      }
	    }
	  }
	}
      }
      if(KFTauCandidates.size()>0){
	KinematicFitTauDecays_.push_back(KFTauCandidates.at(0));
      }
    }
    if (KinematicFitTauDecays_.size() >= minTau_) {
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}


bool KinematicTauProducer::FitKinematicTauCandidate(SelectedKinematicDecay &KFTau,edm::ESHandle<TransientTrackBuilder> &transTrackBuilder_,edm::Handle<reco::GenParticleCollection> &genParticles){
  bool hasasusccessfullfit=false;
  for(unsigned int ambiguity=0; ambiguity<MultiProngTauSolver::NAmbiguity;ambiguity++){
    FitSequencer *kinTauCreator = new ThreeProngTauCreator(transTrackBuilder_, fitParameters_,genParticles);
    int fitStatus = kinTauCreator->create(ambiguity,KFTau);
    if(fitStatus>=1){
      cntSVFound_.at(ambiguity)++;
      ChiSquared chiSquared(kinTauCreator->chi2(0), kinTauCreator->ndf(0));
      if(chiSquared.probability()>0.05)cntSVQC_.at(ambiguity)++;
    }
    
    edm::LogInfo("KinematicTauProducer") <<"KinematicTauProducer::select: fitstatus " << fitStatus ;
    //compute discriminators
    std::map<std::string,bool> discrimValues;
    
    if(fitStatus==2){
      discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFit",true));
      cntLCFit_.at(ambiguity)++;
    }
    else{discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFit",false));}
    bool QC=dicriminatorByKinematicFitQuality(ambiguity,kinTauCreator,fitStatus,KFTau); 
    discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFitQuality",QC));
    if(QC && fitStatus==2){cntLCFitQC_.at(ambiguity)++;;}
    KFTau.SetKinematicFitStatus(ambiguity,discrimValues);
    if(fitStatus==2){
      if(do_BDTTrain_)FillTreeForTraining(ambiguity,kinTauCreator,fitStatus,KFTau);
      saveKinParticles(ambiguity,kinTauCreator,KFTau);
      hasasusccessfullfit=true;
    }
    delete kinTauCreator;
  }
  return hasasusccessfullfit;
}

bool KinematicTauProducer::dicriminatorByKinematicFitQuality(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau){
  //combine a discriminator of loose quality cuts
  //test if fit could create the final decay tree
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality" << std::endl;
  if(fitStatus!=2){ return false;}
  // Configure required paramamters
  reco::PFTau refitPFTau = kinTauCreator->getPFTau();
  std::vector<LorentzVectorParticle> chargedDaughters = kinTauCreator->chargedDaughters();
  std::vector<LorentzVectorParticle> neutralDaughters = kinTauCreator->neutralDaughters();

  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex primaryVtx =KFTau.InitialPrimaryVertexReFit();
  reco::Vertex modifiedPV = KFTau.InitialPrimaryVertexReFitAndRotated();
  reco::Vertex secVtx=ParticleBuilder::GetVertex(kinTauCreator->mother());
  VertexDistance3D vtxdist;
  double vtxSignPVRotSV = vtxdist.distance(modifiedPV, secVtx).significance();
  double vtxSignPVRotPVRed = vtxdist.distance(modifiedPV, primaryVtx).significance();

  // Mass and energy
  double a1Mass = refitPFTau.mass();
  double fraction = refitPFTau.alternatLorentzVect().Et();
  double energyTFraction=-1;
  if(fraction != 0.){energyTFraction = KFTau.PFTauRef()->et()/fraction;}

  // Store quality criteria befora applying
  KFTau.SetQualityCriteria(ambiguity,vtxSignPVRotSV, vtxSignPVRotPVRed, a1Mass, energyTFraction);

  ChiSquared chiSquared(kinTauCreator->chi2(0), kinTauCreator->ndf(0));
  if(chiSquared.probability()<0.05)return false;

  ChiSquared PvtchiSquared(primaryVtx.chi2(),primaryVtx.ndof());
  if(PvtchiSquared.probability()<0.05)return false;
  // Apply selection cuts
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality vtxSignPVRotSV " << vtxSignPVRotSV << std::endl;
  //if ( vtxSignPVRotSV < 2. )return false; // Sig. of secondary vertex
  //    if ( vtxSignPVRotPVRed > 2. )return false; //vertex sig. between modified and initial primary vertex
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality daughters " << chargedDaughters.size() << " " << neutralDaughters.size() << std::endl;
  //if(chargedDaughters.size()!=1 || neutralDaughters.size()!=1) return false; // number of decay products
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality signal cone " << KFTau.PFTauRef()->signalPFChargedHadrCands().size() << std::endl;
  //if(KFTau.PFTauRef()->signalPFChargedHadrCands().size() > 3 ) return false; //tracks in signal cone of initial pftau candidate
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality  Ma1 " << a1Mass << std::endl;
  //if(a1Mass < 0.8) return false; //refitPFTau equals refitted a1 in 3-prong case
  //std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality energyTFraction " << energyTFraction << std::endl;
  //if(energyTFraction == -1.) return false; //energy fraction
  // if(energyTFraction < 0 || energyTFraction > 1.) return false;  //energy fraction
  return true;


}
int KinematicTauProducer::saveKinParticles(unsigned int &ambiguity,FitSequencer * kinTauCreator, SelectedKinematicDecay &KFTau){
  //Get the secondary vertx from fit
  reco::Vertex sec_vertex=ParticleBuilder::GetVertex(kinTauCreator->mother());
  KFTau.SetKFSecondaryVertex(ambiguity,sec_vertex);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Save Kinematic Fit Results and Particles
  //int status = 1;
  std::string name;
  int maxiterations = fitParameters_.getParameter<int>( "maxNbrOfIterations" );
  double mincsum = fitParameters_.getParameter<double>( "maxDelta" );

  SelectedKinematicParticleCollection refitTauDecay;
  kinTauCreator->GetSelectedKinematicParticleList(ambiguity,refitTauDecay);
  KFTau.SetKinematicFitProperties(ambiguity,refitTauDecay, kinTauCreator->Niter(), maxiterations, kinTauCreator->CSum(), mincsum, kinTauCreator->NConstraints(), kinTauCreator->ndf(1), kinTauCreator->chi2(1),1);

  reco::Vertex primaryVtx=KFTau.InitialPrimaryVertexReFit();
  TMatrixT<double> Length(5,1);
  Length(0,0)=sec_vertex.position().x()-primaryVtx.position().x();
  Length(1,0)=sec_vertex.position().y()-primaryVtx.position().y();
  Length(2,0)=sec_vertex.position().z()-primaryVtx.position().z();
  TVector3 l(Length(0,0),Length(1,0),Length(2,0));
  Length(3,0)=l.Phi();
  Length(4,0)=l.Theta();
  TMatrixTSym<double> LengthCov(5);
  for(unsigned int s=0;s<3;s++){
    for(unsigned int t=0;t<3;t++){
      LengthCov(s,t)=primaryVtx.covariance(s,t)+sec_vertex.covariance(s,t);
    }
  }
  TMatrixT<double> Lengthp=MultiProngTauSolver::RotateToTauFrame(Length);
  TMatrixTSym<double> LengthpCov=ErrorMatrixPropagator::PropogateError(&MultiProngTauSolver::RotateToTauFrame,Length,LengthCov);
  double LengthSig=999;
  if(LengthpCov(2,2)!=0)LengthSig=Length(2,0)/LengthpCov(2,2);
  KFTau.SetSecVtxInfo(kinTauCreator->chi2(0), kinTauCreator->ndf(0), Length(2,0), LengthSig);
  return refitTauDecay.size();
}



double KinematicTauProducer::VertexRotationAndSignificance(TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,
                                                                    TVector3 &tauFlghtDirNoCorr,
                                                                    reco::Vertex &pVtx, TLorentzVector &lorentzA1,
                                                                    TVector3 &tauFlghtDir,double &theta0, double &thetaMax){
  TVector3 pv(pVtx.position().x(), pVtx.position().y(), pVtx.position().z());
  TVector3 sv(tmpVtx.position().x(), tmpVtx.position().y(), tmpVtx.position().z());
  tauFlghtDirNoCorr = sv - pv;

  VertexRotation vtxC(lorentzA1);
  thetaMax=fabs(vtxC.calcThetaMax());
  return vtxC.rotatePV(pVtx,tmpVtx,theta0, tauFlghtDir);
}



void KinematicTauProducer::FillTreeForTraining(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau){
  // Configure required paramamters
  reco::PFTau refitPFTau =  kinTauCreator->getPFTau();
  std::vector<LorentzVectorParticle> chargedDaughters = kinTauCreator->chargedDaughters();
  std::vector<LorentzVectorParticle> neutralDaughters = kinTauCreator->neutralDaughters();

  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex primaryVtx=KFTau.InitialPrimaryVertexReFit();
  reco::Vertex modifiedPV=KFTau.InitialPrimaryVertexReFitAndRotated();
  reco::Vertex secVtx=ParticleBuilder::GetVertex(kinTauCreator->mother());
  VertexDistance3D vtxdist;
  double vtxSignPVRotSV = vtxdist.distance(modifiedPV, secVtx).significance();
  double vtxSignPVRotPVRed = vtxdist.distance(modifiedPV, primaryVtx).significance();

  // Mass and energy
  double a1Mass = refitPFTau.mass();
  double fraction = refitPFTau.alternatLorentzVect().Et();
  double energyTFraction=-1;
  if(fraction != 0.){energyTFraction = KFTau.PFTauRef()->et()/fraction;}

  ChiSquared chiSquared(kinTauCreator->chi2(0), kinTauCreator->ndf(0));

  BDT_chiSquared.push_back(TMath::Prob( chiSquared.probability(), 3));
  BDT_energyTFraction.push_back(energyTFraction);
  BDT_vtxSignPVRotSV.push_back(vtxSignPVRotSV);
  BDT_vtxSignPVRotPVRed.push_back(vtxSignPVRotPVRed);
  BDT_a1Mass.push_back(a1Mass);
  BDT_iterations.push_back(kinTauCreator->Niter());
  output_tree->Fill();


}

double KinematicTauProducer::ReturnBDTOutput(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau){
  // Configure required paramamters

  double MvaOut;
  reco::PFTau refitPFTau =  kinTauCreator->getPFTau();
  std::vector<LorentzVectorParticle> chargedDaughters = kinTauCreator->chargedDaughters();
  std::vector<LorentzVectorParticle> neutralDaughters = kinTauCreator->neutralDaughters();

  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex primaryVtx=KFTau.InitialPrimaryVertexReFit();
  reco::Vertex modifiedPV=KFTau.InitialPrimaryVertexReFitAndRotated();
  reco::Vertex secVtx=ParticleBuilder::GetVertex(kinTauCreator->mother());
  VertexDistance3D vtxdist;
  double vtxSignPVRotSV = vtxdist.distance(modifiedPV, secVtx).significance();
  //double vtxSignPVRotPVRed = vtxdist.distance(modifiedPV, primaryVtx).significance();

  // Mass and energy
  double a1Mass = refitPFTau.mass();
  double fraction = refitPFTau.alternatLorentzVect().Et();
  double energyTFraction=-1;
  if(fraction != 0.){energyTFraction = KFTau.PFTauRef()->et()/fraction;}


  ChiSquared chiSquared(kinTauCreator->chi2(0), kinTauCreator->ndf(0));
  if(ambiguity ==0){
    ProbMins3 =TMath::Prob( chiSquared.probability(), 3);
    fracMins =energyTFraction;
    iterMins = kinTauCreator->Niter();
    PVSVMins=vtxSignPVRotSV;
    a1MassMins =a1Mass;
    MvaOut = reader->EvaluateMVA( "BDT_method" );
  }
  else if(ambiguity ==1){
    ProbMins3 =TMath::Prob( chiSquared.probability(), 3);
    fracMins =energyTFraction;
    iterMins = kinTauCreator->Niter();
    PVSVMins=vtxSignPVRotSV;
    a1MassMins =a1Mass;
    MvaOut = reader->EvaluateMVA( "BDT_method" );
  }
  else{//if(ambiguity ==2){
    ProbMins3 =TMath::Prob( chiSquared.probability(), 3);
    fracMins =energyTFraction;
    iterMins = kinTauCreator->Niter();
    PVSVMins=vtxSignPVRotSV;
    a1MassMins =a1Mass;
    MvaOut = reader->EvaluateMVA( "BDT_method" );
  }
  return MvaOut;

}

bool KinematicTauProducer::GetNonTauTracksFromVertex(SelectedKinematicDecay cand,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks){
  const std::vector<reco::TrackRef> tautracks =cand.InitialTrackTriplet();
  const reco::Vertex match=cand.InitialPrimaryVertex();
  // Get track list
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent_->getByLabel(trackCollectionTag_,trackCollection);
  if (!trackCollection.isValid()) {
    edm::LogError("KinematicTauProducer") << "KinematicTauProducer::GetNonTauTracksFromVertex: no track collection found!";
    return false;
  }
  // remove tau tracks and only tracks associated with the vertex
  unsigned int idx = 0;
  for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
    reco::TrackRef tmpRef(trackCollection, idx);
    reco::TrackRef tmpRefForBase=tmpRef;
    if(tmpRef->pt()<17.0){
      bool isTauTrk = false;
      bool fromVertex=false;
      for (std::vector<reco::TrackRef>::const_iterator tauTrk = tautracks.begin(); tauTrk != tautracks.end(); ++tauTrk) {
        if (tmpRef==*tauTrk){isTauTrk = true; break;}
      }
      for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=match.tracks_begin();vtxTrkRef<match.tracks_end();vtxTrkRef++){
        if(match.trackWeight(*vtxTrkRef)>0 ){
          if((*vtxTrkRef)==reco::TrackBaseRef(tmpRefForBase)){fromVertex=true; break;}
        }
      }
      if (!isTauTrk && fromVertex) nonTauTracks.push_back(*iTrk);
    }
  }
  return true;
}



//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
