#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
  fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
  primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
  selectedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "selectedTauCandidates" ) ),//only used to save number of tracks in signal cone of PFTau candidate
  inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) )
{
  produces<reco::PFTauCollection>();
  produces<reco::PFTauDiscriminator>("PFRecoTauDiscriminationByKinematicFit");//boolean per decay whether the fit was successfull
  produces<reco::PFTauDiscriminator>("PFRecoTauDiscriminationByKinematicFitQuality");//boolean per decay whether it passes some major quality cuts
}

KinematicTauProducer::~KinematicTauProducer(){
}

// ------------ method called on each new Event  ------------
bool KinematicTauProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool filterValue = false;
  cnt_++;
  
  std::auto_ptr<reco::PFTauCollection> selected_ = std::auto_ptr<reco::PFTauCollection >(new reco::PFTauCollection);
  reco::PFTauCollection & selected = * selected_;
  
  iEvent_ = &iEvent;
  edm::Handle<reco::VertexCollection> primVtxs;
  iEvent_->getByLabel( primVtx_, primVtxs);
  
  //make copy of initial tau collection with unmodified 4vects
  edm::Handle<reco::PFTauRefVector> usedTaus;
  iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
  selected.insert(selected.begin(), usedTaus->product()->begin(), usedTaus->product()->end());
  
  std::map<int, std::vector<bool> > discrimValues;
  if(primVtxs->size()>=1){
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);
    const reco::Vertex primVtx = primVtxs->front();
    filterValue = select(selected, discrimValues, primVtx);
    if(filterValue) edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Passed";
    else edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::filter Failed"; 
  }
  
  edm::OrphanHandle<reco::PFTauCollection> orphanTaus = iEvent_->put(selected_);
  discriminate(orphanTaus, discrimValues);
	
  if(filterValue) cntFound_++;//found at least 1 refit tau
  
  return filterValue;
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

bool KinematicTauProducer::select(reco::PFTauCollection & selected, std::map<int, std::vector<bool> > & discrimValues, const reco::Vertex & primaryVtx){
  bool success = false;
  discrimValues.clear();
  
  edm::Handle<std::vector<reco::TrackRefVector> > inputCollection;
  iEvent_->getByLabel(inputCollectionTag_, inputCollection);
  edm::Handle<reco::PFTauRefVector> usedTaus;
  iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
  if(inputCollection->size() != usedTaus->size()){
    edm::LogError("KinematicTauProducer")<<"KinematicTauProducer::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<selectedTauCandidatesTag_.label()<<"("<<usedTaus->size()<<")";
    return false;
  }
  edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::select: Size of usedTaus: "<< usedTaus->size() << " Size of Track Collection: " << inputCollection->size()  ;
  unsigned int index = 0;
  //TransientTrackBuilder trkBuilder = *transTrackBuilder_;
  KinematicTauCreator *kinTauCrtr = new ThreeProngTauCreator(transTrackBuilder_, fitParameters_);
  for(std::vector<reco::TrackRefVector>::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks, ++index) {
    std::vector<reco::TrackRef> input;
    for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
    int fitStatus = kinTauCrtr->create(primaryVtx, input);
    edm::LogInfo("KinematicTauProducer")<<"KinematicTauProducer::select: fitstatus " << fitStatus << " Tracks " << input.size() 
					<< " Vertex (" << primaryVtx.position().x() << "," << primaryVtx.position().y() << "," << primaryVtx.position().z() << ")";       
    //save std::map<int, std::vector<bool> >, stores association between position in tauCollection and one bool for every discriminator
    reco::PFTauRef tauRef = usedTaus->at(index);
    discrimValues.insert(std::make_pair(tauRef.index(), std::vector<bool>()));
    discrimValues.find(tauRef.index())->second.push_back(fitStatus);
    discrimValues.find(tauRef.index())->second.push_back(dicriminatorByKinematicFitQuality(kinTauCrtr, fitStatus, tauRef, primaryVtx));
    
    if(fitStatus==1){
      success = true;
      //modify tau in selected list
      reco::PFTau refitPFTau = kinTauCrtr->getPFTau();//this is only the visible part of the tau momentum!
      selected.at(tauRef.index()).setP4(refitPFTau.p4());//this is only the visible part of the tau momentum!
      selected.at(tauRef.index()).setalternatLorentzVect(kinTauCrtr->getKinematicTau().p4());//this is the refitted part of the tau momentum including the neutrino!
      selected.at(tauRef.index()).setVertex(refitPFTau.vertex()); //this is the rotated primary vertex
    }
  }
  
  delete kinTauCrtr;
  
  return success; //at least one tau was fitted
}
bool KinematicTauProducer::dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr, const int & fitStatus, const reco::PFTauRef & tauRef, const reco::Vertex & primaryVtx){
  //combine a discriminator of loose quality cuts
  //	bool debug = true;
  
  //test if fit could create the final decay tree
  if(!fitStatus){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad fit status! cntFound is "<<cntFound_;
    return false;
  }
  reco::PFTau refitPFTau = kinTauCrtr->getPFTau();//this is only the visible part of the refitted tau momentum!
  refitPFTau.setalternatLorentzVect(kinTauCrtr->getKinematicTau().p4());//this is the whole refitted tau momentum including the neutrino!
  std::vector<math::XYZTLorentzVector> chargedDaughters = kinTauCrtr->getRefittedChargedDaughters();
  std::vector<math::XYZTLorentzVector> neutralDaughters = kinTauCrtr->getRefittedNeutralDaughters();
  
  //chi2prob
  ChiSquared chiSquared(kinTauCrtr->chi2(), kinTauCrtr->ndf());
  //	if( fabs(TMath::Prob(kinTauCrtr->chi2(), kinTauCrtr->ndf()) - chiSquared.probability()) > 0.00001) printf("KinematicTauProducer::dicriminatorByKinematicFitQuality: tested probs differ. TMath %f, CMSSW %f\n", TMath::Prob(kinTauCrtr->chi2(), kinTauCrtr->ndf()), chiSquared.probability());
  if( chiSquared.probability() < 0.03 ){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad chi2prob of "<<chiSquared.probability()<<"! cntFound is "<<cntFound_;
    return false;
  }
  
  //vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
  reco::Vertex modifiedPV = kinTauCrtr->getModifiedPrimaryVertex();
  VertexState secVtx(kinTauCrtr->getKinematicTree()->currentDecayVertex()->position(), kinTauCrtr->getKinematicTree()->currentDecayVertex()->error());
  VertexDistance3D vtxdist;
  if ( vtxdist.distance(modifiedPV, secVtx).significance() < 2. ){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad SV sign of "<<vtxdist.distance(modifiedPV, secVtx).significance()<<"! cntFound is "<<cntFound_;
    return false;
  }
  //vertex significance between modified and initial primary vertex
  if ( vtxdist.distance(modifiedPV, primaryVtx).significance() > 2. ){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad PV rot of "<<vtxdist.distance(modifiedPV, primaryVtx).significance()<<"! cntFound is "<<cntFound_;
    return false;
  }
  
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!
  if(chargedDaughters.size()!=3 || neutralDaughters.size()!=1){
    edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality:WARNING!!! KinematicTauProducer assumes a tau decay into three pions and neutrino but recieved "<<chargedDaughters.size()<<" charged and "<<neutralDaughters.size()<<" neutral daughters!";
    return false;
  }
  
  //tracks in signal cone of initial pftau candidate
  if ( tauRef->signalPFChargedHadrCands().size() > 3 ){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad track size of "<<tauRef->signalPFChargedHadrCands().size()<<"! cntFound is "<<cntFound_;
    return false;
  }
  
  //a1 mass
  if( refitPFTau.mass() < 0.8){
    //if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad a1 mass of "<<refitPFTau.mass()<<"! cntFound is "<<cntFound_;
    return false;//refitPFTau equals refitted a1 in 3-prong case
  }
  
  //energy fraction
  double fraction = refitPFTau.alternatLorentzVect().Et();
  if( fraction == 0.){
    edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality:WARNING!!! Bad energy in alternatLorentzVect of 0! visible energy is "<<refitPFTau.et()<<".";		
    return false;
  }	
  fraction = tauRef->et()/fraction;
  if(fraction < 0 || fraction > 1.){
    //		if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: Bad energy ratio of "<<fraction<<"! cntFound is "<<cntFound_;
    return false;
  }
  
  //	if(debug) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality: counter is "<<cntFound_<<" values are: chi2prob "<<chiSquared.probability()<<", SV sep "<<vtxdist.distance(modifiedPV, secVtx).significance()<<", PV rot"<<vtxdist.distance(modifiedPV, primaryVtx).significance()<<", trks size "<<tauRef->signalPFChargedHadrCands().size()<<", a1 mass "<<refitPFTau.mass()<<", fraction"<<fraction;
  
  return true;
}
void KinematicTauProducer::discriminate(const edm::OrphanHandle<reco::PFTauCollection> & collection, const std::map<int, std::vector<bool> > & discrimValues){
  std::auto_ptr<reco::PFTauDiscriminator> discrKinFit = std::auto_ptr<reco::PFTauDiscriminator>(new reco::PFTauDiscriminator(reco::PFTauRefProd(collection)));
  std::auto_ptr<reco::PFTauDiscriminator> discrKinFitQual = std::auto_ptr<reco::PFTauDiscriminator>(new reco::PFTauDiscriminator(reco::PFTauRefProd(collection)));
  
  for(std::map<int, std::vector<bool> >::const_iterator iter = discrimValues.begin(); iter != discrimValues.end(); ++iter){
    reco::PFTauRef tauRef(collection, iter->first);
    if(iter->second.size() < 2){
      edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::discriminate: ERROR! Bad discriminator size. This tau will be skipped.";
      continue;
    }
    discrKinFit->setValue(iter->first, iter->second.at(0));
    discrKinFitQual->setValue(iter->first, iter->second.at(1));
  }
  
  iEvent_->put(discrKinFit, "PFRecoTauDiscriminationByKinematicFit");
  iEvent_->put(discrKinFitQual, "PFRecoTauDiscriminationByKinematicFitQuality");
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
