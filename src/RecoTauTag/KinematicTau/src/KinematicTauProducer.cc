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
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag"))
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

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

  filterValue = select(KinematicFitTauDecays_,daughterCollection);
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
}

void KinematicTauProducer::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("KinematicTau")<<"--> [KinematicTauProducer] asks for >= 1 kinTau per event. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection &KinematicFitTauDecays_,reco::RecoChargedCandidateCollection & daughterCollection){
  bool success = false;

  edm::Handle<SelectedKinematicDecayCollection > KinematicTauCandidates;
  iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidates);
  
  KinematicTauCreator *kinTauCreator = new ThreeProngTauCreator(transTrackBuilder_, fitParameters_);
  for(unsigned int i=0;i<KinematicTauCandidates->size();i++){
    SelectedKinematicDecay KFTau=KinematicTauCandidates->at(i);
    int fitStatus = kinTauCreator->create(KFTau);
    edm::LogInfo("KinematicTauProducer") <<"KinematicTauProducer::select: fitstatus " << fitStatus ;
    
    //compute discriminators
    std::map<std::string,bool> discrimValues;
    discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFit",fitStatus));
    discrimValues.insert(std::pair<std::string,bool>("PFRecoTauDiscriminationByKinematicFitQuality",dicriminatorByKinematicFitQuality(kinTauCreator,fitStatus,KFTau)));

    if(fitStatus==1){
      //modify tau in selected list
      KinematicFitTauDecays_.push_back(KFTau);
      std::vector<reco::TrackRef> usedTracks = kinTauCreator->getSelectedTracks();
      saveSelectedTracks(usedTracks, daughterCollection);
      saveKinParticles(kinTauCreator,KinematicFitTauDecays_,discrimValues);
      success =true;
    }
  }
  return success; //at least one tau was fitted
}

bool KinematicTauProducer::dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr, const int & fitStatus, SelectedKinematicDecay &KFTau){
  //combine a discriminator of loose quality cuts
  //test if fit could create the final decay tree
  if(!fitStatus)return false;
  // Configure required paramamters
  reco::PFTau refitPFTau = kinTauCrtr->getPFTau();//this is only the visible part of the refitted tau momentum!
  refitPFTau.setalternatLorentzVect(kinTauCrtr->getKinematicTau().p4());//this is the whole refitted tau momentum including the neutrino!
  std::vector<math::XYZTLorentzVector> chargedDaughters = kinTauCrtr->getRefittedChargedDaughters();
  std::vector<math::XYZTLorentzVector> neutralDaughters = kinTauCrtr->getRefittedNeutralDaughters();
  //chi2prob
  ChiSquared chiSquared(kinTauCrtr->chi2(), kinTauCrtr->ndf());
  if( chiSquared.probability() < 0.03 )return false;
  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex primaryVtx =KFTau.PrimaryVertexReFit();
  reco::PFTauRef tauRef=KFTau.PFTauRef();
  reco::Vertex modifiedPV = kinTauCrtr->getModifiedPrimaryVertex();
  VertexState secVtx(kinTauCrtr->getKinematicTree()->currentDecayVertex()->position(), kinTauCrtr->getKinematicTree()->currentDecayVertex()->error());
  VertexDistance3D vtxdist;
  double fraction = refitPFTau.alternatLorentzVect().Et();
  fraction = tauRef->et()/fraction;
  
  // Apply selection cuts
  if ( vtxdist.distance(modifiedPV, secVtx).significance() < 2. )return false; // Sig. of secondary vertex
  if ( vtxdist.distance(modifiedPV, primaryVtx).significance() > 2. )return false; //vertex sig. between modified and initial primary vertex
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!
  if(chargedDaughters.size()!=3 || neutralDaughters.size()!=1) return false; // number of decay products
  if(tauRef->signalPFChargedHadrCands().size() > 3 ) return false; //tracks in signal cone of initial pftau candidate
  if(refitPFTau.mass() < 0.8) return false; //refitPFTau equals refitted a1 in 3-prong case
  if(fraction == 0.) return false; //energy fraction
  if(fraction < 0 || fraction > 1.) return false;  //energy fraction
  return true;
}


int KinematicTauProducer::saveKinParticles(const KinematicTauCreator * kinTauCrtr, SelectedKinematicDecayCollection &KinFitTau, std::map<std::string, bool> tauDiscriminators){
  RefCountedKinematicTree tree = kinTauCrtr->getKinematicTree();
  KinematicConstrainedVertexFitter *kcvFitter = kinTauCrtr->getFitter();
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
  int ambiguityCnt = -1;
  int maxiterations = fitParameters_.getParameter<int>( "maxNbrOfIterations" );
  double mincsum = fitParameters_.getParameter<double>( "maxDelta" );
  reco::RecoChargedCandidateRef emptyCandRef;//pions will be filled later on by correctReferences()      
  if(tree->currentParticle()->currentState().particleCharge() != 0){
    name = std::string("tau");
  }
  else{
    LogTrace("KinematicTauProducer")<<"KinematicTauProducer::saveKinParticles: neutral tau detected. tau skipped.";
    return 0;
  }
  reco::PFTauRef tauRef=KinFitTau.back().PFTauRef(); // Fix this is not the real inital state
  SelectedKinematicParticleCollection refitTauDecay;
  refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), status, name, ambiguityCnt, emptyCandRef) );
  refitTauDecay.back().setInitialState(TLorentzVector(tauRef->px(), tauRef->py(), tauRef->pz(), tauRef->energy()), kinTauCrtr->getModifiedPrimaryVertex());//initial tau state consists of rotated primVtx (+init. err) and the pftau par. 
  int constraints = tree->currentParticle()->degreesOfFreedom();
  std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
  for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter) {
    if((*iter)->currentState().particleCharge() != 0){
      name = std::string("pion");
    }else{
      name = std::string("neutrino");
    }
    refitTauDecay.push_back( SelectedKinematicParticle(*iter, status, name, ambiguityCnt, emptyCandRef) );
  }
  
 if(refitTauDecay.size() != 5) LogTrace("KinematicTauProducer")<<"KinematicTauProducer::saveSelectedTracks:Saved only "<<refitTauDecay.size()<<" refitted particles.";
 else{
   KinFitTau.back().SetKinematicFitProperties(refitTauDecay, iterations, maxiterations, csum, mincsum, constraints, kinTauCrtr->ndf(), kinTauCrtr->chi2(),tauDiscriminators);
   setMissingQualityCriteria(KinFitTau, kinTauCrtr);
 }
 return refitTauDecay.size();
}



void KinematicTauProducer::saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection){
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

void KinematicTauProducer::correctReferences(SelectedKinematicDecayCollection & KinFitTaus, const edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands){
  unsigned index = 0;
  std::vector<reco::RecoChargedCandidateRef> newRefs;
  for(reco::RecoChargedCandidateCollection::const_iterator iter = orphanCands->begin(); iter != orphanCands->end(); ++iter, ++index){
    reco::RecoChargedCandidateRef ref(orphanCands, index);
    newRefs.push_back(ref);
  }
  index = 0;
  for(SelectedKinematicDecayCollection::iterator decay = KinFitTaus.begin(); decay != KinFitTaus.end(); ++decay){
    std::vector< SelectedKinematicParticle* > daughters;
    decay->modifiableChargedDaughters(daughters);
    for(std::vector<SelectedKinematicParticle*>::iterator particle = daughters.begin(); particle != daughters.end(); ++particle){
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

void KinematicTauProducer::setMissingQualityCriteria(SelectedKinematicDecayCollection &KinFitTau, const KinematicTauCreator * kinTauCrtr){
  // Get saved Values 
  reco::PFTauRef tauRef=KinFitTau.back().PFTauRef();
  reco::Vertex   primaryVtx=KinFitTau.back().PrimaryVertexReFit();
  reco::PFTau    refitPFTau = kinTauCrtr->getPFTau();//this is only the visible part of the refitted tau momentum!
  refitPFTau.setalternatLorentzVect(kinTauCrtr->getKinematicTau().p4());//this is the whole refitted tau momentum including the neutrino!
  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex modifiedPV = kinTauCrtr->getModifiedPrimaryVertex();
  VertexState secVtx(kinTauCrtr->getKinematicTree()->currentDecayVertex()->position(), kinTauCrtr->getKinematicTree()->currentDecayVertex()->error());
  VertexDistance3D vtxdist;
  double vtxSignPVRotSV = vtxdist.distance(modifiedPV, secVtx).significance();
  //vertex significance between modified and initial primary vertex  
  double vtxSignPVRotPVRed = vtxdist.distance(modifiedPV, primaryVtx).significance();
  //a1 mass
  double a1Mass = refitPFTau.mass();
  //energy fraction  
  double energyTFraction = refitPFTau.alternatLorentzVect().Et();
  if( energyTFraction == 0.){
    edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality:WARNING!!! Bad energy in alternatLorentzVect of 0! visible energy is "<<refitPFTau.et()<<".";
    energyTFraction = -1.;
  }
  else {
    energyTFraction = tauRef->et()/energyTFraction;
  }
  KinFitTau.back().setMissingQualityCriteria(vtxSignPVRotSV, vtxSignPVRotPVRed, a1Mass, energyTFraction);
}


//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
