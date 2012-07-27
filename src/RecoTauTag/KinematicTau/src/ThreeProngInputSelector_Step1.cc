#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step1.h"

ThreeProngInputSelector_Step1::ThreeProngInputSelector_Step1(const edm::ParameterSet& iConfig):
  KinematicTauTools(),
  inputCollectionTag_(iConfig.getParameter<edm::InputTag>("tauCandidates")),
  trackCollectionTag_(iConfig.getParameter<edm::InputTag>("tracks")),
  minTau_(iConfig.getUntrackedParameter<unsigned int>("minTau", 1)) //filter returns true if more/equal than minTau_ taus were selected
{
  iConfig_ = iConfig;
  produces<int>("flag"); //0=invalid, 1=valid
  produces<reco::TrackCollection>("NonTauTracks"); //save collection of tracks not belonging to any tau candidate
  produces<vVVTrackRef>("ThreeProngCombinations"); //save collection of vVVTrackRef for each tau candidate
}


ThreeProngInputSelector_Step1::~ThreeProngInputSelector_Step1() {
}

// ------------ method called on each new Event  ------------
bool ThreeProngInputSelector_Step1::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  cnt_++;
  iEvent_ = &iEvent;
  
  std::auto_ptr<reco::TrackCollection> pNonTauTracks = std::auto_ptr<reco::TrackCollection >(new reco::TrackCollection);
  reco::TrackCollection & nonTauTracks = * pNonTauTracks;
  
  std::auto_ptr<vVVTrackRef > pThreeProngCombis = std::auto_ptr<vVVTrackRef >(new vVVTrackRef);
  vVVTrackRef & threeProngCombis = * pThreeProngCombis;

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  Set_TransientTrackBuilder(transTrackBuilder);
  
  bool filterValue = select(nonTauTracks, threeProngCombis);
  
  iEvent_->put(pNonTauTracks, "NonTauTracks");
  iEvent_->put(pThreeProngCombis, "ThreeProngCombinations");
  
  std::auto_ptr<int> pFlag = std::auto_ptr<int>(new int);
  int & flag = * pFlag;
  if (filterValue) {
    flag = 1;
  } 
  else {
    flag = 0;
  }
  iEvent_->put(pFlag, "flag");
  
  return filterValue;
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeProngInputSelector_Step1::beginJob() {
  cnt_ = 0;
  cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeProngInputSelector_Step1::endJob() {
  float ratio = 0.0;
  if (cnt_ != 0) ratio = (float)cntFound_/cnt_;
  edm::LogVerbatim("ThreeProngInputSelector_Step1") << "--> [ThreeProngInputSelector_Step1] found at least " << minTau_ << " 3-prongs per event. Efficiency: " << cntFound_ << "/" << cnt_ << " = " << std::setprecision(4) << ratio*100.0 << "%";
}

bool ThreeProngInputSelector_Step1::select(reco::TrackCollection & nonTauTracks, vVVTrackRef & threeProngCombis) {
  std::vector<reco::TrackRef> tautracks;
  
  //load input collection from InputTrackSelector
  //each tau candidate stores its signal cone tracks
  edm::Handle<std::vector<reco::TrackRefVector> > inputCollection;
  iEvent_->getByLabel(inputCollectionTag_, inputCollection);
  
  if (!inputCollection.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: no std::vector<reco::TrackRefVector> found!";
    return false;
  }
  
  for (std::vector<reco::TrackRefVector>::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks) {
    std::vector<reco::TrackRef> input;
    for (reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
    threeProngCombis.push_back(choose3Prongs(input)); //make combinations of exact 3 tracks, choose3Prongs can return an empty vector which will be deleted later on (together with its tauRef)
  }
  
  //save all signal cone tracks of all candidates ONCE into tautracks without duplicates
  //keep duplicates in duplicateTracks and delete combi if all three tracks are also in ONE other combi (of another candidate)
  std::vector<std::vector<std::vector<reco::TrackRef> > >::iterator candidates, passedCandidates;
  std::vector<std::vector<reco::TrackRef> > ::iterator triplets, passedTriplets;
  std::vector<reco::TrackRef>::const_iterator tracks, trackrefs;
  for (candidates = threeProngCombis.begin(); candidates != threeProngCombis.end(); ++candidates) {
    for (triplets = candidates->begin(); triplets != candidates->end();) {
      bool erasedTriplet = false;
      std::vector<reco::TrackRef> duplicateTracks;
      for (tracks = triplets->begin(); tracks != triplets->end(); ++tracks) {
	bool exists = false;
	for (trackrefs = tautracks.begin(); trackrefs != tautracks.end(); ++trackrefs) {
	  if (*tracks == *trackrefs) {
	    exists = true;
	    break;
	  }
	}
	if (!exists) {
	  tautracks.push_back(*tracks);
	} 
	else {
	  duplicateTracks.push_back(*tracks);   
	}
      }
      if (duplicateTracks.size() == 3) erasedTriplet = removeDuplicateTriplets(duplicateTracks, threeProngCombis, candidates, triplets);
      if (!erasedTriplet) ++triplets;
    }
  }
  
  //load general track collection and substract tau tracks from it
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent_->getByLabel(trackCollectionTag_, trackCollection);
  
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
    
  if (threeProngCombis.size() >= minTau_) {
    LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: " << threeProngCombis.size() << " tau candidate(s) selected.";
    cntFound_++;
    return true;
  } 
  else {
    LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: Warning! Only " << threeProngCombis.size() << " tau candidate(s) reconstructed. Skip Event.";
    return false;
  }
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step1);
