#include "RecoTauTag/KinematicTau/interface/InputTrackSelector.h"

InputTrackSelector::InputTrackSelector(const edm::ParameterSet& iConfig):
  tauType_( iConfig.getUntrackedParameter<std::string>("tauType", "hps") ),
  trkCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauDaughterTracks" ) ),  // Restrict tau daughters to origin from a certain track collection
  //  vtxtrackCollectionTag_(iConfig.getParameter<edm::InputTag>("vtxtracks")),            // Track collection for the vertices
  primVtx_(iConfig.getParameter<edm::InputTag>("primVtx")),
  minTracks_( iConfig.getParameter<unsigned int>("minTracks") ),                    // Only tau candidates with more/equal than minTracks are selected
  minTau_( iConfig.getUntrackedParameter<unsigned int>("minTau", 1) ),              // Filter returns true if more/equal than minTau_ taus were selected
  nTauPerVtx_(iConfig.getUntrackedParameter<unsigned int>("nTauPerVtx", 0)),         // Number of taus per vertex: 0=no requirement; 1 = 1 tau per vertex only; 2 = 2 taus per vertex only
  minTauPt_( iConfig.getUntrackedParameter<double>("minTauPt", 0.) ),               // Ignore pftaus below this pt threshold
  TauVtxList_( iConfig.getUntrackedParameter< std::vector<std::string> >("NonTauTracks") )
{
  produces<std::vector<std::vector<SelectedKinematicDecay> > >("PreKinematicDecaysStep1");
  for(unsigned int i=0; i<TauVtxList_.size(); i++){
    produces<reco::TrackCollection>(TauVtxList_.at(i)); //save collection of tracks not belonging to any tau candidate
  }
}


InputTrackSelector::~InputTrackSelector(){
}

// ------------ method called on each new Event  ------------
void InputTrackSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  cnt_++;
  iEvent_ = &iEvent;
  std::auto_ptr<std::vector<std::vector<SelectedKinematicDecay> > > KFCandidates_ = std::auto_ptr<std::vector<std::vector<SelectedKinematicDecay> > >(new std::vector<std::vector<SelectedKinematicDecay> >);
  std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates = *KFCandidates_;

  std::vector<reco::TrackCollection> NonTauTracksLists_;
  for(unsigned int i=0;i<TauVtxList_.size(); i++){
    NonTauTracksLists_.push_back(reco::TrackCollection());
  }

  select(KFCandidates,NonTauTracksLists_);

  iEvent_->put(KFCandidates_,"PreKinematicDecaysStep1");
  for(unsigned int i=0;i<TauVtxList_.size(); i++){
    std::auto_ptr<reco::TrackCollection > NonTauTracks_ = std::auto_ptr<reco::TrackCollection >(new reco::TrackCollection);
    *NonTauTracks_ =NonTauTracksLists_.at(i);
    iEvent_->put(NonTauTracks_,TauVtxList_.at(i));
  }

}

// ------------ method called once each job just before starting event loop  ------------
void InputTrackSelector::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void InputTrackSelector::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  edm::LogVerbatim("InputTrackSelector")<<"--> [InputTrackSelector] found at least "<<minTau_<<" tau candidate(s) of type "<<tauType_<<" (with at least "<<minTracks_<<" tracks and pt > "<<minTauPt_<<") per event. Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

bool InputTrackSelector::select(std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates,std::vector<reco::TrackCollection> &NonTauTracksLists_){
  bool found = false;
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtx_,primaryVertexCollection);
  
  if(primaryVertexCollection->size()>0){
    edm::Handle<reco::PFTauCollection> inputCollection;
    iEvent_->getByLabel(tauType_+"PFTauProducer", inputCollection);
    
    edm::Handle<reco::TrackCollection> trkCollection;
    iEvent_->getByLabel(trkCollectionTag_, trkCollection);
    trkCollectionID_ = trkCollection.id();
    unsigned int p=0;
    for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
      reco::PFTauRef thePFTau(inputCollection, iPFTau);
      // filter candidates
      if (thePFTau->pt() < minTauPt_) continue;
      reco::TrackRefVector tauDaughters = getPFTauDaughters(thePFTau);
      if(tauDaughters.size()>=minTracks_){
	// Add KF Taus Here
	unsigned int tau_idx=KFCandidates.size();
	KFCandidates.push_back(std::vector<SelectedKinematicDecay>());
	std::vector<reco::TrackRef> input;
	for (reco::TrackRefVector::iterator trk =tauDaughters.begin(); trk!=tauDaughters.end(); ++trk) input.push_back(*trk);
	std::vector<std::vector<reco::TrackRef> > combis=choose3Prongs(input);
	for(unsigned int i=0;i<combis.size();i++){
	  if(NonTauTracksLists_.size()==1 && NonTauTracksLists_.size()==TauVtxList_.size() && nTauPerVtx_==0){
	    KFCandidates.at(tau_idx).push_back(SelectedKinematicDecay(SelectedKinematicDecay::ThreePion,thePFTau,combis.at(i),primaryVertexCollection->front(),TauVtxList_.at(0),nTauPerVtx_));
	  }
	  else if(NonTauTracksLists_.size()>1 && TauVtxList_.size()==NonTauTracksLists_.size() && nTauPerVtx_==1){
	    //Code to find best vertex would go here
	    KFCandidates.at(tau_idx).push_back(SelectedKinematicDecay(SelectedKinematicDecay::ThreePion,thePFTau,combis.at(i),primaryVertexCollection->front(),TauVtxList_.at(p),nTauPerVtx_));
	    p++;
	  }
	  else{
	    edm::LogWarning("InputTrackSelector")<<"InputTrackSelector::select WARNING!!! Number of allowed verticies out of range size: " <<  TauVtxList_.size() << " value: " << p
						   <<"Skipping Tau candidate ";
	  }
	}
      }
    }
      
    // remove duplicates
    std::vector<reco::TrackRef> tautracks;
    for(unsigned int i=0;i<KFCandidates.size();i++){
      for(unsigned int j=0;j<KFCandidates.at(i).size();j++){
	std::vector<reco::TrackRef> Triplet_ij=KFCandidates.at(i).at(j).InitalTrackTriplet();
	std::vector<reco::TrackRef>::iterator track_ij;
	for (track_ij = Triplet_ij.begin(); track_ij != Triplet_ij.end(); track_ij++) tautracks.push_back(*track_ij);
	for(unsigned int k=0;k<KFCandidates.size();k++){
	  for(unsigned int l=0;l<KFCandidates.at(k).size();l++){
	    if(!(i==k && j==l)){
	      unsigned int duplicates=0;
	      std::vector<reco::TrackRef> Triplet_kl=KFCandidates.at(k).at(l).InitalTrackTriplet();
	      std::vector<reco::TrackRef>::iterator track_kl;
	      for (track_ij = Triplet_ij.begin(); track_ij != Triplet_ij.end(); track_ij++){
		for (track_kl = Triplet_kl.begin(); track_kl != Triplet_kl.end(); track_kl++){
		  if (*track_ij==*track_kl) {
		    duplicates++;
		  }
		}
	      }
	      if(duplicates==3){KFCandidates.at(k).erase(KFCandidates.at(k).begin()+l); l--;}
	    }
	  }
	}
      }
    }
    
    // Get Vertex Tracks List
    if(NonTauTracksLists_.size()==1 && nTauPerVtx_==0){
      GetNonTauTracks(iEvent_,trkCollectionTag_,NonTauTracksLists_.at(0),tautracks);
    }
    else if(nTauPerVtx_==1){
      unsigned int p=0;
      for(unsigned int i=0;i<KFCandidates.size();i++){
	for(unsigned int j=0;j<KFCandidates.at(i).size();j++){
	  if(p<NonTauTracksLists_.size()){
	    GetNonTauTracksFromVertex(KFCandidates.at(i).at(j),trkCollectionTag_,NonTauTracksLists_.at(p));
	  }
	  p++;
	}
      }
    }
  }
  unsigned int ntaus(0);
  for(unsigned int i=0;i<KFCandidates.size();i++){
    for(unsigned int j=0;j<KFCandidates.at(i).size();j++){
      ntaus++;
    }
  }
  if(ntaus>= minTau_){
    cntFound_++;
    found = true;
    LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: "<< ntaus <<" tau candidate(s) reconstructed.";
  }
  else LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: Only "<< ntaus <<" tau candidate(s) reconstructed. Skip Evt.";
    
  return found;
}

reco::TrackRefVector InputTrackSelector::getPFTauDaughters(reco::PFTauRef &PFTau){
  reco::TrackRefVector trkVct;
  const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands(); //candidates in signal cone 
  //isolationPFChargedHadrCands stores tracks in isolation/veto cone
  for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
    if(iter->get()->trackRef().isNonnull()) {
      if ((*iter)->trackRef().id() != trkCollectionID_) {
	// ignore tracks that do not origin from the desired track collection (e.g. ignore conversionStepTracks) for now
	const edm::Provenance & prov = iEvent_->getProvenance((*iter)->trackRef().id());
	edm::LogInfo("InputTrackSelector")<<"InputTrackSelector::getPFTauDaughters: Skip PFTau daughter from collection "<<prov.moduleLabel();
	continue;
      }
      trkVct.push_back( (*iter)->trackRef() );
    }
  }
  return trkVct;
}

bool  InputTrackSelector::GetNonTauTracks(edm::Event *iEvent_,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks, std::vector<reco::TrackRef> &tautracks){
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent_->getByLabel(trackCollectionTag_,trackCollection);
  if (!trackCollection.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: no track collection found!";
    return false;
  }
  unsigned int idx = 0;
  for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
    reco::TrackRef tmpRef(trackCollection, idx);
    bool isTauTrk = false;
    for (std::vector<reco::TrackRef>::const_iterator tauTrk = tautracks.begin(); tauTrk != tautracks.end(); ++tauTrk) {
      if (tmpRef == *tauTrk) {
        isTauTrk = true;
        break;
      }
    }
    if (!isTauTrk) nonTauTracks.push_back(*iTrk);
  }
  return true;
}


bool InputTrackSelector::GetNonTauTracksFromVertex(SelectedKinematicDecay cand,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks){
  const std::vector<reco::TrackRef> tautracks =cand.InitalTrackTriplet();
  const reco::Vertex match=cand.InitalPrimaryVertex();
  // Get track list
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent_->getByLabel(trackCollectionTag_,trackCollection);
  if (!trackCollection.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: no track collection found!";
    return false;
  }
  // remove tau tracks and only tracks associated with the vertex
  unsigned int idx = 0;
  for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
    reco::TrackRef tmpRef(trackCollection, idx);
    reco::TrackRef tmpRefForBase=tmpRef;
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
  return true;
}



// Compute the charge of the triplet
bool InputTrackSelector::sumCharge(const std::vector<reco::TrackRef> & input){
  return (1== abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge()));
}

// compute all permutations
std::vector<std::vector<reco::TrackRef> > InputTrackSelector::permuteCombinations(const std::vector<reco::TrackRef> & vect){
  std::vector<std::vector<reco::TrackRef> > combis;
  std::vector<reco::TrackRef>::const_iterator iter1, iter2, iter3;
  for (iter1 = vect.begin(); iter1 != vect.end(); ++iter1) {
    iter2 = iter1;
    ++iter2;
    for (; iter2 != vect.end(); ++iter2) {
      iter3 = iter2;
      ++iter3;
      for (; iter3 != vect.end(); ++iter3) {
	std::vector<reco::TrackRef> newCombi;
        newCombi.push_back(*iter1);
        newCombi.push_back(*iter2);
        newCombi.push_back(*iter3);
        combis.push_back(newCombi);
      }
    }
  }
  return combis;
}

// Select 3 prong based on mass and charge       
std::vector<std::vector<reco::TrackRef> > InputTrackSelector::choose3Prongs(std::vector<reco::TrackRef> &input){
  sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);                                                                                      
  std::vector<std::vector<reco::TrackRef> > combis = permuteCombinations(input);
  for (std::vector<std::vector<reco::TrackRef> >::iterator iter = combis.begin(); iter != combis.end();) {  
    if (!sumCharge(*iter)) { 
      iter = combis.erase(iter); 
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong charge sum. " << combis.size() << " combis left.";
      continue;  
    }                              
    double massA1 = getInvariantMass(*iter);  
    if (massA1 > 2.0 || massA1 < 3*PMH.Get_piMass()) { //soft upper value     
      iter = combis.erase(iter);
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong mass. " << combis.size() << " combis left.";  
      continue;
    } 
    ++iter; 
  }  
  return combis; 
}                   

double InputTrackSelector::getInvariantMass(std::vector<reco::TrackRef> &tracks){
  double SumPx(0),SumPy(0),SumPz(0),SumE(0); 
  for(unsigned int i=0; i<tracks.size(); i++){
    SumPx += tracks.at(i)->px();              
    SumPy += tracks.at(i)->py();              
    SumPz += tracks.at(i)->pz();              
    SumE += sqrt(pow(tracks.at(i)->p(),2)+pow(PMH.Get_piMass(),2)); 
  }
  return sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2));
}

template <typename T>  bool InputTrackSelector::cmpPt(const T & a, const T & b){
  return a->pt() > b->pt();
}


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);




