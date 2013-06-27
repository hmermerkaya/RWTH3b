#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step2.h"

ThreeProngInputSelector_Step2::ThreeProngInputSelector_Step2(const edm::ParameterSet & iConfig):
  KinematicTauTools(),
  threeProngCollectionTag_(iConfig.getParameter<edm::InputTag>("threeProngs")),
  selectedTauCandidatesTag_(iConfig.getParameter<edm::InputTag>("selectedTauCandidates")),
  primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
  minTau_(iConfig.getUntrackedParameter<unsigned int>("minTau", 1)), //filter returns true if more/equal than minTau_ taus were selected
  minVtxTracks_(iConfig.getUntrackedParameter("minVtxTracks", int(3))),
  maxChi2ndf_(iConfig.getUntrackedParameter("maxChi2ndf", double(10.0)))
{
    iConfig_ = iConfig;
    produces<int>("flag"); //0=invalid, 1=valid
    produces<std::vector<reco::TrackRefVector> >("InputTracks"); //save collection of vector<reco::CandidateRef> for each tau cand
    produces<reco::PFTauRefVector>("InputTauRefs"); //needed to fill in unfit KinematicParticle later on
    produces<reco::VertexCollection>("primVtx"); //has to be vector. save one of length one
}


ThreeProngInputSelector_Step2::~ThreeProngInputSelector_Step2() {
}

// ------------ method called on each new Event  ------------
bool ThreeProngInputSelector_Step2::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  cnt_++;
  iEvent_ = &iEvent;

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  Set_TransientTrackBuilder(transTrackBuilder);
  
  std::auto_ptr<std::vector<reco::TrackRefVector> > pSelected = std::auto_ptr<std::vector<reco::TrackRefVector> >(new std::vector<reco::TrackRefVector>);
  std::vector<reco::TrackRefVector> & selected = * pSelected;
  
  std::auto_ptr<reco::PFTauRefVector> pPFTauRef = std::auto_ptr<reco::PFTauRefVector>(new reco::PFTauRefVector);
  reco::PFTauRefVector & PFTauRefs = * pPFTauRef;
  
  std::auto_ptr<reco::VertexCollection> pPrimaryVertex = std::auto_ptr<reco::VertexCollection>(new reco::VertexCollection);
  reco::VertexCollection & primaryVertex = * pPrimaryVertex;
  
  bool filterValue = select(selected, PFTauRefs, primaryVertex);
  
  iEvent_->put(pSelected, "InputTracks");
  iEvent_->put(pPFTauRef, "InputTauRefs");
  iEvent_->put(pPrimaryVertex, "primVtx");
  
  std::auto_ptr<int> pFlag = std::auto_ptr<int>(new int);
  int & flag = * pFlag;
  if (filterValue) flag = 1;
  else flag = 0;
  iEvent_->put(pFlag, "flag");
  
  return filterValue;
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeProngInputSelector_Step2::beginJob() {
  cnt_ = 0;
  cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeProngInputSelector_Step2::endJob() {
  float ratio = 0.0;
  if (cnt_ != 0) ratio = (float)cntFound_/cnt_;
  edm::LogVerbatim("ThreeProngInputSelector_Step2") << "--> [ThreeProngInputSelector_Step2] found at least " << minTau_ << " 3-prongs per event. Efficiency: " << cntFound_ << "/" << cnt_ << " = " << std::setprecision(4) << ratio*100.0 << "%";
}

bool ThreeProngInputSelector_Step2::select(std::vector<reco::TrackRefVector> & selected, reco::PFTauRefVector & taurefs, reco::VertexCollection & selectedPrimaryVertex) {
  //load three-prong collection from step 1
  edm::Handle<vVVTrackRef > threeProngCombinations;
  iEvent_->getByLabel(threeProngCollectionTag_, threeProngCombinations);
  
  if (!threeProngCombinations.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: no three-prong collection found!";
    return false;
  }
  
  vVVTrackRef threeProngCombis = *threeProngCombinations;
  
  //load pftau references
  edm::Handle<reco::PFTauRefVector> inputTauCollection;
  iEvent_->getByLabel(selectedTauCandidatesTag_, inputTauCollection);
  taurefs = *inputTauCollection;
  
  if (!inputTauCollection.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: no inputTauCollection found!";
    return false;
  }
  if (threeProngCombinations->size() != inputTauCollection->size()) {
    edm::LogError("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Bad input collections. Size mismatch between " << threeProngCollectionTag_.label() << "(" << threeProngCombinations->size() << ") and " << selectedTauCandidatesTag_.label() << "(" << inputTauCollection->size() << ")";
    return false;
  }
  
  //load primary vertex collection created w/o tau tracks
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtxTag_, primaryVertexCollection);
  
  if (!primaryVertexCollection.isValid()) {
    edm::LogError("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: no primary-vertex collection found!";
    return false;
  }
  
  reco::VertexCollection primaryVertices = *primaryVertexCollection;
  
  //add additional primary-vertex cleaning here
  if (primaryVertices.size() > 0) {
    selectedPrimaryVertex.push_back(primaryVertices.front());
  }
  
  if (selectedPrimaryVertex.size() > 0) {
    reco::PFTauRefVector::iterator pfiter = taurefs.begin();
    for (vVVTrackRef::iterator candidates = threeProngCombis.begin(); candidates != threeProngCombis.end(); ++candidates) {
      bool success = choose3bestTracks(selected, *candidates, selectedPrimaryVertex.front());
      if (!success) {
	pfiter = taurefs.erase(pfiter);
      } 
      else {
	++pfiter;
      }
    }
    if (selected.size() >= minTau_) {
      LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: " << selected.size() << " tau candidate(s) reconstructed.";
      cntFound_++;
      return true;
    } 
    else {
      LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Warning! Only " << selected.size() << " tau candidate(s) reconstructed. Skip Event.";
      return false;
    }
  }
  LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Unable to calculate a new primary vertex. No tau candidates reconstructed.";
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step2);
