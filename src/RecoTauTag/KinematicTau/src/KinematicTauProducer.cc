#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"
#include "RecoTauTag/KinematicTau/interface/ParticleBuilder.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
  fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
  primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag")),
  VertexTags_(iConfig.getUntrackedParameter< std::vector<std::string> >("VertexTags")),
  TauVtxList_(iConfig.getUntrackedParameter< std::vector<std::string> >("NonTauTracks")),
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
  if(filterValue) cntFound_++;//found at least 1 refit tau
  iEvent_->put(KinematicFitTauDecays,"KinematicFitTau");
}

void KinematicTauProducer::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;
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
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] asks for >= 1 kinTau per event. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
  
  if(do_BDTTrain_){
    output->Write();
    output->Close();
  }

}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,const edm::EventSetup& iSetup){
  std::cout << "KinematicTauProducer::select" << std::endl;
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
    std::cout << "found vertex" << std::endl;
    edm::Handle<std::vector<std::vector<SelectedKinematicDecay> > > KinematicTauCandidate;
    iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidate);
    std::cout << "N candidates List " << KinematicTauCandidate->size() << std::endl;
    for(unsigned int i=0; i<KinematicTauCandidate->size();i++){
      std::vector<SelectedKinematicDecay>  KFTauCandidates;
      std::cout << "Kinematic Fit " << i << std::endl;
      std::cout << "N candidates" << KinematicTauCandidate->at(i).size() << std::endl;
      for(unsigned int j=0; j<KinematicTauCandidate->at(i).size();j++){
	std::cout << "Kinematic Fit " << i << " " << j << std::endl;
	SelectedKinematicDecay KTau=KinematicTauCandidate->at(i).at(j);
        SecondaryVertexHelper SVH(transTrackBuilder_,KTau);
        if(SVH.hasSecondaryVertex()){
	  std::cout << "Kinematic Fit " << i << " " << j << " has SV" << std::endl;
          TString vertexName=KTau.PrimaryVertexReFitCollectionTag();
          TString VTag;
          for(unsigned int v=0;v<VertexTags_.size() && VertexTags_.size()==TauVtxList_.size();v++){
            if(vertexName==TauVtxList_.at(v)) VTag=VertexTags_.at(v);
          }
	  edm::Handle<reco::VertexCollection > CurrentTauPrimaryVtx;
          iEvent_->getByLabel(edm::InputTag(VTag.Data()),CurrentTauPrimaryVtx);
	  std::cout << "Kinematic Fit check PV " << std::endl;
          if(!CurrentTauPrimaryVtx.isValid()) continue;
          if(CurrentTauPrimaryVtx->size()==0) continue;
	  std::cout << "Kinematic Fit has PV " << std::endl;
	  reco::Vertex primaryVertexReFit=CurrentTauPrimaryVtx->front();
	  reco::Vertex primaryVertexReFitAndRotated=primaryVertexReFit;
          TVector3 tauFlghtDirNoCorr;
          TVector3 tauFlghtDir;
          TLorentzVector a1_p4=SVH.Initial_a1_p4();
          double initThetaGJ,ThetaMax;
          TransientVertex SecondaryVertex=SVH.InitialSecondaryVertex();
	  std::vector<reco::TransientTrack> RefittedTracks=SVH.InitialRefittedTracks();
          double s = VertexRotationAndSignificance(SecondaryVertex,RefittedTracks,tauFlghtDirNoCorr,primaryVertexReFitAndRotated,a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	  if(/*(s<sigcut_  || fabs(initThetaGJ)<fabs(ThetaMax)) &&*/ s>=0 ){
	    if(fabs(tauFlghtDirNoCorr.Eta())<etacut_ || fabs(tauFlghtDir.Eta())<etacut_ ){
	      KTau.SetInitialVertexProperties(primaryVertexReFit,primaryVertexReFitAndRotated,SVH.InitialRefittedTracks(),SVH.InitialSecondaryVertex());
	      KTau.SetInitialKinematics(tauFlghtDirNoCorr,SVH.Initial_pions(),a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	      std::cout << "before Fit "<< std::endl;
	      bool FitOK=FitKinematicTauCandidate(KTau,transTrackBuilder_,genParticles);
	      std::cout << "after Fit "<< std::endl;
	      if(FitOK){
		std::cout << "good Fit "<< std::endl;
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
      LogTrace("KinematicTauProducer") << "KinematicTauProducer::select: " << KinematicFitTauDecays_.size() << " tau candidate(s) reconstructed.";
      cntFound_++;
      return true;
    }
    else {
      LogTrace("KinematicTauProducer") << "KinematicTauProducer::select: Warning! Only " << KinematicFitTauDecays_.size() << " tau candidate(s) reconstructed. Skip Event.";
      return false;
    }
  }
  LogTrace("KinematicTauProducer") << "KinematicTauProducer::select: Unable to calculate a new primary vertex. No tau candidates reconstructed.";
  return false;
}


bool KinematicTauProducer::FitKinematicTauCandidate(SelectedKinematicDecay &KFTau,edm::ESHandle<TransientTrackBuilder> &transTrackBuilder_,edm::Handle<reco::GenParticleCollection> &genParticles){
  bool hasasusccessfullfit=false;
  for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){
    FitSequencer *kinTauCreator = new ThreeProngTauCreator(transTrackBuilder_, fitParameters_,genParticles);
    int fitStatus = kinTauCreator->create(ambiguity,KFTau);
    edm::LogInfo("KinematicTauProducer") <<"KinematicTauProducer::select: fitstatus " << fitStatus ;
    //compute discriminators
    std::map<std::string,bool> discrimValues;
    
    discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFit",fitStatus));
    discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFitQuality",dicriminatorByKinematicFitQuality(ambiguity,kinTauCreator,fitStatus,KFTau)));
    KFTau.SetKinematicFitStatus(ambiguity,discrimValues);
    std::cout << "Fit status " << fitStatus << std::endl;
    if(fitStatus==1){
      std::cout << "Fit status A " << std::endl;
      if(do_BDTTrain_)FillTreeForTraining(ambiguity,kinTauCreator,fitStatus,KFTau);
      std::cout << "Fit status B " << std::endl;
      std::cout<<"----------- BDT  "<<ReturnBDTOutput(ambiguity,kinTauCreator,fitStatus,KFTau) <<std::endl;
      saveKinParticles(ambiguity,kinTauCreator,KFTau);
      std::cout << "Fit status C " << std::endl;
      hasasusccessfullfit=true;
    }
    delete kinTauCreator;
  }
  return hasasusccessfullfit;
}

bool KinematicTauProducer::dicriminatorByKinematicFitQuality(unsigned int &ambiguity,FitSequencer *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau){
  //combine a discriminator of loose quality cuts
  //test if fit could create the final decay tree
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality" << std::endl;
  if(!fitStatus){std::cout << "Fit Failed" << std::endl; return false;}
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

  ChiSquared chiSquared(kinTauCreator->chi2(), kinTauCreator->ndf());

  //if( chiSquared.probability() < 0.03 )return false;
  // Apply selection cuts
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality vtxSignPVRotSV " << vtxSignPVRotSV << std::endl;
  //if ( vtxSignPVRotSV < 2. )return false; // Sig. of secondary vertex
  //    if ( vtxSignPVRotPVRed > 2. )return false; //vertex sig. between modified and initial primary vertex
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality daughters " << chargedDaughters.size() << " " << neutralDaughters.size() << std::endl;
  //if(chargedDaughters.size()!=1 || neutralDaughters.size()!=1) return false; // number of decay products
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality signal cone " << KFTau.PFTauRef()->signalPFChargedHadrCands().size() << std::endl;
  //if(KFTau.PFTauRef()->signalPFChargedHadrCands().size() > 3 ) return false; //tracks in signal cone of initial pftau candidate
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality  Ma1 " << a1Mass << std::endl;
  //if(a1Mass < 0.8) return false; //refitPFTau equals refitted a1 in 3-prong case
  std::cout << "KinematicTauProducer::dicriminatorByKinematicFitQuality energyTFraction " << energyTFraction << std::endl;
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
  KFTau.SetKinematicFitProperties(ambiguity,refitTauDecay, kinTauCreator->Niter(), maxiterations, kinTauCreator->CSum(), mincsum, kinTauCreator->NConstraints(), kinTauCreator->ndf(), kinTauCreator->chi2(),1);
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

  ChiSquared chiSquared(kinTauCreator->chi2(), kinTauCreator->ndf());

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


  ChiSquared chiSquared(kinTauCreator->chi2(), kinTauCreator->ndf());
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

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
