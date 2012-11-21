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



KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
  fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" ))
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
  std::auto_ptr<reco::RecoChargedCandidateCollection> daughterCollection_ = std::auto_ptr<reco::RecoChargedCandidateCollection>(new reco::RecoChargedCandidateCollection);
  reco::RecoChargedCandidateCollection & daughterCollection = * daughterCollection_;

  filterValue = select(KinematicFitTauDecays_,daughterCollection,iSetup);
  if(filterValue) edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Passed";
  else edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Failed";

  edm::OrphanHandle<reco::RecoChargedCandidateCollection> orphanCands = iEvent_->put(daughterCollection_,"KinematicFitTauDaughters");
  correctReferences(KinematicFitTauDecays_, orphanCands) ;//has to be called before put(selected_,"SelectedKinematicParticles")!!!
  iEvent_->put(KinematicFitTauDecays,"KinematicFitTau");

  if(filterValue) cntFound_++;//found at least 1 refit tau
}

void KinematicTauProducer::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;

   output = new TFile("ForTrain.root","RECREATE");
   output_tree = new TTree("t","t");
   output_tree->Branch("BDT_vtxSignPVRotSV",&BDT_vtxSignPVRotSV);
   output_tree->Branch("BDT_vtxSignPVRotPVRed",&BDT_vtxSignPVRotPVRed);
   output_tree->Branch("BDT_a1Mass",&BDT_a1Mass);
   output_tree->Branch("BDT_energyTFraction",&BDT_energyTFraction);
   output_tree->Branch("BDT_chiSquared",&BDT_chiSquared);

}

void KinematicTauProducer::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] asks for >= 1 kinTau per event. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
  output->Write();
  output->Close();

}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,reco::RecoChargedCandidateCollection & daughterCollection,const edm::EventSetup& iSetup){
  bool success = false;

  edm::Handle<SelectedKinematicDecayCollection > KinematicTauCandidates;
  iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidates);

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent_->getByLabel(gensrc_, genParticles);

  for(unsigned int i=0;i<KinematicTauCandidates->size();i++){

    SelectedKinematicDecay KFTau=KinematicTauCandidates->at(i);
    bool hasasusccessfullfit=false;

    for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity;ambiguity++){

      KinematicTauCreator *kinTauCreator = new ThreeProngTauCreator(transTrackBuilder_, fitParameters_,genParticles);

      int fitStatus = kinTauCreator->create(ambiguity,KFTau);

      edm::LogInfo("KinematicTauProducer") <<"KinematicTauProducer::select: fitstatus " << fitStatus ;
      if(fitStatus==1)kinTauCreator->getRefittedChargedDaughters();
      //compute discriminators
      std::map<std::string,bool> discrimValues;

      discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFit",fitStatus));
      discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFitQuality",dicriminatorByKinematicFitQuality(ambiguity,kinTauCreator,fitStatus,KFTau)));
      KFTau.SetKinematicFitStatus(ambiguity,discrimValues);
      std::cout<< "fitStatus "<< fitStatus << "by QC  " << dicriminatorByKinematicFitQuality(ambiguity,kinTauCreator,fitStatus,KFTau) <<std::endl;
      if(fitStatus==1){
 
	saveKinParticles(ambiguity,kinTauCreator,KFTau);
	saveSelectedTracks(kinTauCreator->getSelectedTracks(),daughterCollection);
	hasasusccessfullfit=true;
	success =true;
      }
      delete kinTauCreator;
    }
    if(hasasusccessfullfit)KinematicFitTauDecays_.push_back(KFTau);
  }
  return success; //at least one tau was fitted
}

bool KinematicTauProducer::dicriminatorByKinematicFitQuality(unsigned int &ambiguity,const KinematicTauCreator *kinTauCreator, const int & fitStatus, SelectedKinematicDecay &KFTau){
  //combine a discriminator of loose quality cuts
  //test if fit could create the final decay tree
  if(!fitStatus)return false;
  // Configure required paramamters
  reco::PFTau refitPFTau = kinTauCreator->getPFTau();//this is only the visible part of the refitted tau momentum!
  refitPFTau.setalternatLorentzVect(kinTauCreator->getKinematicTau().p4());//this is the whole refitted tau momentum including the neutrino!
  std::vector<math::XYZTLorentzVector> chargedDaughters = kinTauCreator->getRefittedChargedDaughters();
  std::vector<math::XYZTLorentzVector> neutralDaughters = kinTauCreator->getRefittedNeutralDaughters();

  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex primaryVtx =KFTau.InitialPrimaryVertexReFit();
  reco::Vertex modifiedPV = kinTauCreator->getModifiedPrimaryVertex();
  VertexState secVtx(kinTauCreator->getKinematicTree()->currentDecayVertex()->position(), kinTauCreator->getKinematicTree()->currentDecayVertex()->error());
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
  BDT_chiSquared =chiSquared.probability();
  BDT_energyTFraction =energyTFraction;
  BDT_vtxSignPVRotSV =vtxSignPVRotSV;
  BDT_vtxSignPVRotPVRed =vtxSignPVRotPVRed;
  BDT_a1Mass =a1Mass;
  output_tree->Fill();

  //if( chiSquared.probability() < 0.03 )return false;
  // Apply selection cuts
    if ( vtxSignPVRotSV < 2. )return false; // Sig. of secondary vertex
    //    if ( vtxSignPVRotPVRed > 2. )return false; //vertex sig. between modified and initial primary vertex
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!
  if(chargedDaughters.size()!=1 || neutralDaughters.size()!=1) return false; // number of decay products
  if(KFTau.PFTauRef()->signalPFChargedHadrCands().size() > 3 ) return false; //tracks in signal cone of initial pftau candidate
  if(a1Mass < 0.8) return false; //refitPFTau equals refitted a1 in 3-prong case
  if(energyTFraction == -1.) return false; //energy fraction
  if(energyTFraction < 0 || energyTFraction > 1.) return false;  //energy fraction
  return true;


}


int KinematicTauProducer::saveKinParticles(unsigned int &ambiguity,const KinematicTauCreator * kinTauCreator, SelectedKinematicDecay &KFTau){
  //Get the secondary vertx from fit
  kinTauCreator->getKinematicTree()->movePointerToTheTop();
  reco::Vertex sec_vertex;
  RefCountedKinematicParticle the_top = kinTauCreator->getKinematicTree()->currentParticle();
  if (the_top->currentState().isValid()){
    KinematicVertex the_vertex = (*kinTauCreator->getKinematicTree()->currentDecayVertex());
    if(the_vertex.vertexIsValid() ){
      sec_vertex = reco::Vertex(reco::Vertex::Point(the_vertex.vertexState().position()),the_vertex.vertexState().error().matrix_new(),
				the_vertex.chiSquared(), the_vertex.degreesOfFreedom(),0);
    }
  }
  KFTau.SetKFSecondaryVertex(ambiguity,sec_vertex);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Save Kinematic Fit Results and Particles
  RefCountedKinematicTree tree = kinTauCreator->getKinematicTree();
  NumericalKinematicConstrainedFitter *kcvFitter = kinTauCreator->getFitter();
  try{
    tree->movePointerToTheTop();
  }
  catch(VertexException){
    edm::LogWarning("KinematicTauProducer")<<"KinematicTree::movePointerToTheTop; tree is empty! -- Event skipped.";
    return false;
  }
  int iterations = kcvFitter->getNit();
  float csum = kcvFitter->getCSum();
  
  int status = 1;
  std::string name;
  int maxiterations = fitParameters_.getParameter<int>( "maxNbrOfIterations" );
  double mincsum = fitParameters_.getParameter<double>( "maxDelta" );
  reco::RecoChargedCandidateRef emptyCandRef;//pions will be filled later on by correctReferences()      
  if(tree->currentParticle()->currentState().particleCharge() != 0){
    name = std::string("tau");std::cout<<"tau"<<std::endl;
  }
  else{
    LogTrace("KinematicTauProducer")<<"KinematicTauProducer::saveKinParticles: neutral tau detected. tau skipped.";
    return 0;
  }
  reco::PFTauRef tauRef=KFTau.PFTauRef(); // Fix this is not the real initial state
  SelectedKinematicParticleCollection refitTauDecay;
  refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), status, name, ambiguity, emptyCandRef) );
  refitTauDecay.back().setInitialState(TLorentzVector(tauRef->px(), tauRef->py(), tauRef->pz(), tauRef->energy()), kinTauCreator->getModifiedPrimaryVertex());//initial tau state consists of rotated primVtx (+init. err) and the pftau par. 
  int constraints = tree->currentParticle()->degreesOfFreedom();
  std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();

  for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter) {

    if((*iter)->currentState().particleCharge() != 0){
      name = std::string("a1");std::cout<<"a1"<<std::endl;
    }else{
      name = std::string("neutrino");std::cout<<"neutrino"<<std::endl;
    }
  
    refitTauDecay.push_back( SelectedKinematicParticle(*iter, status, name, ambiguity, emptyCandRef) );
  }


  std::vector<RefCountedKinematicParticle> Pions = kinTauCreator->getPions();
  for (std::vector<RefCountedKinematicParticle>::iterator itr=Pions.begin(); itr!=Pions.end(); ++itr) {  
    if((*itr)->currentState().particleCharge() != 0){
      name = std::string("pion");std::cout<<"pion"<<std::endl;
    }
     refitTauDecay.push_back( SelectedKinematicParticle(*itr, status, name, ambiguity, emptyCandRef) );
  }

  if(refitTauDecay.size() != 6) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::saveKinParticles Invalid number of SelectedKinematicParticles saveSelectedTracks:Saved only "<<refitTauDecay.size()<<" refitted particles.";
 else{
   KFTau.SetKinematicFitProperties(ambiguity,refitTauDecay, iterations, maxiterations, csum, mincsum, constraints, kinTauCreator->ndf(), kinTauCreator->chi2());
 }
 return refitTauDecay.size();
}



void KinematicTauProducer::saveSelectedTracks(const std::vector<reco::TrackRef> usedTracks, reco::RecoChargedCandidateCollection & daughterCollection){
  //set up TrackToCandidate converter              
  edm::ParameterSet pioncfg;
  pioncfg.addParameter("particleType", 211);
  converter::TrackToCandidate tk2cand(pioncfg);
  for (std::vector<reco::TrackRef>::const_iterator trk = usedTracks.begin(); trk != usedTracks.end(); ++trk) {
    reco::RecoChargedCandidate tmpCand;
    tk2cand.convert(*trk, tmpCand); //convert tracks to candidates                                                                                                                                                                           
    daughterCollection.push_back(tmpCand);
  }
}


void KinematicTauProducer::correctReferences(SelectedKinematicDecayCollection & KFTaus, const edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands){
  unsigned index = 0;
  std::vector<reco::RecoChargedCandidateRef> newRefs;
  for(reco::RecoChargedCandidateCollection::const_iterator iter = orphanCands->begin(); iter != orphanCands->end(); ++iter, ++index){
    reco::RecoChargedCandidateRef ref(orphanCands, index);
    newRefs.push_back(ref);
  }
  index = 0;
  for(SelectedKinematicDecayCollection::iterator decay = KFTaus.begin(); decay != KFTaus.end(); ++decay){
    std::vector< SelectedKinematicParticle* > daughters;

    decay->modifiableChargedDaughters(daughters);
    for(std::vector<SelectedKinematicParticle*>::iterator particle = daughters.begin(); particle != daughters.end(); ++particle){
      std::cout << "Solution: " << (*particle)->p4().M() << " " <<  (*particle)->p4().Px() << " " <<  (*particle)->p4().Py() << " " <<  (*particle)->p4().Pz() << std::endl; 
      if(index>=newRefs.size()){
	edm::LogError("KinematicTauProducer")<<"evt "<<iEvent_->id().event()<<" KinematicTauProducer::correctReferences: Bad selection size! index="<<index<<", refs="<<newRefs.size();
        throw 111;
      }
      if(*particle!=NULL){
        (*particle)->setCandRef(newRefs.at(index));
      }else edm::LogError("KinematicTauProducer")<<"KinematicTauProducer::correctReferences: Reference not modified!!!(at index="<<index<<")";
      index++;
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
