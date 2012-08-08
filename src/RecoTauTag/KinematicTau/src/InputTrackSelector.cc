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

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  Set_TransientTrackBuilder(transTrackBuilder);

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
	vVTrackRef combis=choose3Prongs(input);
	if(NonTauTracksLists_.size()==1 || nTauPerVtx_==0){
	  for(unsigned int i=0;i<combis.size();i++){
	    KFCandidates.at(tau_idx).push_back(SelectedKinematicDecay(SelectedKinematicDecay::ThreePion,thePFTau,combis.at(i),primaryVertexCollection->front(),TauVtxList_.at(0),nTauPerVtx_));
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
    if(NonTauTracksLists_.size()==1)GetNonTauTracks(iEvent_,trkCollectionTag_,NonTauTracksLists_.at(0),tautracks);
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
    //LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::getPFTauDaughters: PFTau daughter pt "<<iter->get()->pt()<<", eta "<<iter->get()->eta()<<", vtx("<<iter->get()->vx()<<","<<iter->get()->vy()<<","<<iter->get()->vz()<<")";
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


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);
