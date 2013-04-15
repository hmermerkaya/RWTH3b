#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step1.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

ThreeProngInputSelector_Step1::ThreeProngInputSelector_Step1(const edm::ParameterSet& iConfig):
  tauType_( iConfig.getUntrackedParameter<std::string>("tauType", "hps") ),
  trkCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauDaughterTracks" ) ),  // Restrict tau daughters to origin from a certain track collection
  //  vtxtrackCollectionTag_(iConfig.getParameter<edm::InputTag>("vtxtracks")),            // Track collection for the vertices
  primVtx_(iConfig.getParameter<edm::InputTag>("primVtx")),
  minTracks_( iConfig.getParameter<unsigned int>("minTracks") ),                    // Only tau candidates with more/equal than minTracks are selected
  minTau_( iConfig.getUntrackedParameter<unsigned int>("minTau", 1) ),              // Filter returns true if more/equal than minTau_ taus were selected
  minTauPt_( iConfig.getUntrackedParameter<double>("minTauPt", 10.) ),               // Ignore pftaus below this pt threshold
  TauEtaCut_( iConfig.getUntrackedParameter<double>("TauEtaCut", 2.0) )
{
  produces<std::vector<std::vector<SelectedKinematicDecay> > >("PreKinematicDecaysStep1");
}


ThreeProngInputSelector_Step1::~ThreeProngInputSelector_Step1(){
}

// ------------ method called on each new Event  ------------
void ThreeProngInputSelector_Step1::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  cnt_++;
  iEvent_ = &iEvent;

  std::auto_ptr<std::vector<std::vector<SelectedKinematicDecay> > > KFCandidates_ = std::auto_ptr<std::vector<std::vector<SelectedKinematicDecay> > >(new std::vector<std::vector<SelectedKinematicDecay> >);
  std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates = *KFCandidates_;

  select(KFCandidates);

  iEvent_->put(KFCandidates_,"PreKinematicDecaysStep1");
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeProngInputSelector_Step1::beginJob(){
  cnt_ = 0;
  cntFound_ = 0;
  n3track=0;
  n3prong=0; 
  nhps3prong=0;
  npftaus=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeProngInputSelector_Step1::endJob(){
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  std::cout <<"--> [ThreeProngInputSelector_Step1] found at least "<<minTau_<<" tau candidate(s) of type "<<tauType_<<" (with at least "<<minTracks_<<" tracks and pt > "<<minTauPt_<<") per event. Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%" << std::endl;

  std::cout <<"[ThreeProngInputSelector_Step1] found at least 3 Prong Candidates  Efficiency: " << n3prong <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)n3prong)/((float)cnt_)*100.0<<"%" << std::endl;
  std::cout <<"[ThreeProngInputSelector_Step1] found at least tracks  Efficiency: " << n3track <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)n3track)/((float)cnt_)*100.0<<"%"<< std::endl;
  std::cout <<"[ThreeProngInputSelector_Step1] found PFT with HPS 3 prongs  Efficiency: " << nhps3prong <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)nhps3prong)/((float)cnt_)*100.0<<"%"<< std::endl;
  std::cout <<"[ThreeProngInputSelector_Step1] found at least 1 PFT Efficiency: " << npftaus <<"/"<<cnt_<<" = "<<std::setprecision(4)<< ((float)npftaus)/((float)cnt_)*100.0<<"%"<< std::endl;
}

bool ThreeProngInputSelector_Step1::select(std::vector<std::vector<SelectedKinematicDecay> > &KFCandidates){
  bool found = false;
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtx_,primaryVertexCollection);
  if(primaryVertexCollection->size()>0){
    edm::Handle<reco::PFTauCollection> inputCollection;
    iEvent_->getByLabel(tauType_+"PFTauProducer", inputCollection);
    
    edm::Handle<reco::TrackCollection> trkCollection;
    iEvent_->getByLabel(trkCollectionTag_, trkCollection);
    trkCollectionID_ = trkCollection.id();
    bool has3tracks=false;
    bool hasHPS3prong=false;
    bool has3prongcandidate=false;
    for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
     reco::PFTauRef thePFTau(inputCollection, iPFTau);
      // filter candidates
      if(thePFTau->pt() < minTauPt_) continue;
      if(fabs(thePFTau->eta())>TauEtaCut_)continue;
      reco::TrackRefVector tauDaughters = getPFTauDaughters(thePFTau);
      if(tauDaughters.size()>=3) has3tracks=true;
      if(thePFTau->decayMode()==10)hasHPS3prong=true;
      if(tauDaughters.size()>=minTracks_){
	// Add KF Taus Here
	unsigned int tau_idx=KFCandidates.size();
	KFCandidates.push_back(std::vector<SelectedKinematicDecay>());
	std::vector<reco::TrackRef> input;
	for (reco::TrackRefVector::iterator trk =tauDaughters.begin(); trk!=tauDaughters.end(); ++trk){/*if((*trk)->algo()==reco::TrackBase::iter0 || (*trk)->algo()==reco::TrackBase::iter3)*/ input.push_back(*trk);}
	std::vector<std::vector<reco::TrackRef> > combis=choose3Prongs(input);
	for(unsigned int i=0;i<combis.size();i++){
	  has3prongcandidate=true;
	  KFCandidates.at(tau_idx).push_back(SelectedKinematicDecay(SelectedKinematicDecay::ThreePion,thePFTau,combis.at(i),primaryVertexCollection->front(),"notused",1));
	}
      }
    }
    if(has3tracks)n3track++;
    if(hasHPS3prong)nhps3prong++;
    if(has3prongcandidate)n3prong++;
    if(inputCollection->size()>0)npftaus++;
    // remove duplicates
    std::vector<reco::TrackRef> tautracks;
    for(unsigned int i=0;i<KFCandidates.size();i++){
      for(unsigned int j=0;j<KFCandidates.at(i).size();j++){
	std::vector<reco::TrackRef> Triplet_ij=KFCandidates.at(i).at(j).InitialTrackTriplet();
	std::vector<reco::TrackRef>::iterator track_ij;
	for (track_ij = Triplet_ij.begin(); track_ij != Triplet_ij.end(); track_ij++) tautracks.push_back(*track_ij);
	for(unsigned int k=0;k<KFCandidates.size();k++){
	  for(unsigned int l=0;l<KFCandidates.at(k).size();l++){
	    if(!(i==k && j==l)){
	      unsigned int duplicates=0;
	      std::vector<reco::TrackRef> Triplet_kl=KFCandidates.at(k).at(l).InitialTrackTriplet();
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
    LogTrace("ThreeProngInputSelector_Step1")<<"evt "<<iEvent_->id().event()<<" ThreeProngInputSelector_Step1::select: "<< ntaus <<" tau candidate(s) reconstructed.";
  }
  else LogTrace("ThreeProngInputSelector_Step1")<<"evt "<<iEvent_->id().event()<<" ThreeProngInputSelector_Step1::select: Only "<< ntaus <<" tau candidate(s) reconstructed. Skip Evt.";
  return found;
}

reco::TrackRefVector ThreeProngInputSelector_Step1::getPFTauDaughters(reco::PFTauRef &PFTau){
  reco::TrackRefVector trkVct;
  const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands(); //candidates in signal cone 
  if(cands.size()==0) return trkVct;
  unsigned int j=0;
  double ptmax=0;
  for(unsigned int i=0;i<cands.size();i++){
    if(cands.at(i)->trackRef().isNonnull()) {
      trkVct.push_back(cands.at(i)->trackRef());
    }
  }


  /*
  for(unsigned int i=0;i<cands.size();i++){
    if(cands.at(i)->trackRef().isNonnull()) {
      if(cands.at(i)->trackRef()->pt()>ptmax){ptmax=cands.at(i)->trackRef()->pt(); j=i;}
    }
  }
  
  if(!cands.at(j)->trackRef().isNonnull() || cands.at(j)->trackRef().id() != trkCollectionID_) return trkVct;
  trkVct.push_back( cands.at(j)->trackRef() );
  TLorentzVector pi;
  pi.SetXYZM(trkVct.at(0)->px(),trkVct.at(0)->py(),trkVct.at(0)->pz(),PDGInfo::pi_mass());
  const reco::PFCandidateRefVector &  outercands = PFTau->isolationPFChargedHadrCands();
  for(unsigned int i=0;i<cands.size();i++){
    if(i!=j && cands.at(i)->trackRef().isNonnull() && cands.at(i)->trackRef().id() == trkCollectionID_){
      TLorentzVector LV;
      LV.SetXYZM(cands.at(i)->trackRef()->px(),cands.at(i)->trackRef()->py(),cands.at(i)->trackRef()->pz(),PDGInfo::pi_mass());
      LV+=pi;
      if(LV.M()<PDGInfo::tau_mass())trkVct.push_back(cands.at(i)->trackRef());
    }
  }
  for(unsigned int i=0;i<outercands.size();i++){
    if(outercands.at(i)->trackRef().isNonnull() && outercands.at(i)->trackRef().id() == trkCollectionID_){
      TLorentzVector LV; 
      LV.SetXYZM(outercands.at(i)->trackRef()->px(),outercands.at(i)->trackRef()->py(),outercands.at(i)->trackRef()->pz(),PDGInfo::pi_mass());
      LV+=pi;
      if(LV.M()<PDGInfo::tau_mass())trkVct.push_back(outercands.at(i)->trackRef());
    }
  }
  */
  return trkVct;
}

bool  ThreeProngInputSelector_Step1::GetNonTauTracks(edm::Event *iEvent_,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks, std::vector<reco::TrackRef> &tautracks){
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


bool ThreeProngInputSelector_Step1::GetNonTauTracksFromVertex(SelectedKinematicDecay cand,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks){
  const std::vector<reco::TrackRef> tautracks =cand.InitialTrackTriplet();
  const reco::Vertex match=cand.InitialPrimaryVertex();
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



// Compute the charge of the triplet
bool ThreeProngInputSelector_Step1::sumCharge(const std::vector<reco::TrackRef> & input){
  return (1== abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge()));
}

// compute all permutations
std::vector<std::vector<reco::TrackRef> > ThreeProngInputSelector_Step1::permuteCombinations(const std::vector<reco::TrackRef> & vect){
  std::vector<std::vector<reco::TrackRef> > combis;
  std::vector<reco::TrackRef>::const_iterator iter1, iter2, iter3;
  if(vect.begin()!=vect.end()){
    iter1 = vect.begin(); 
    iter2 = iter1;
    iter2++;
    for (; iter2 != vect.end(); iter2++) {
      iter3 = iter2;
      iter3++;
      for (; iter3 != vect.end(); iter3++) {
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
std::vector<std::vector<reco::TrackRef> > ThreeProngInputSelector_Step1::choose3Prongs(std::vector<reco::TrackRef> &input){
  sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);                                                                                      
  std::vector<std::vector<reco::TrackRef> > combis = permuteCombinations(input);
  for (std::vector<std::vector<reco::TrackRef> >::iterator iter = combis.begin(); iter != combis.end();) {  
    if (!sumCharge(*iter)) { 
      iter = combis.erase(iter); 
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong charge sum. " << combis.size() << " combis left.";
      continue;  
    }                              
    double massA1 = getInvariantMass(*iter);  
    double pt = getSumPt(*iter);
    if (massA1 > 2.0 || massA1 < 3*PMH.Get_piMass() || pt<12.0) { //soft upper value     
      iter = combis.erase(iter);
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong mass. " << combis.size() << " combis left.";  
      continue;
    } 
    ++iter; 
  }  
  return combis; 
}                   

double ThreeProngInputSelector_Step1::getSumPt(std::vector<reco::TrackRef> &tracks){
  double SumPx(0),SumPy(0);
  for(unsigned int i=0; i<tracks.size(); i++){
    SumPx += tracks.at(i)->px();
    SumPy += tracks.at(i)->py();
  }
  return sqrt(pow(SumPx,2)+pow(SumPy,2));
}

double ThreeProngInputSelector_Step1::getInvariantMass(std::vector<reco::TrackRef> &tracks){
  double SumPx(0),SumPy(0),SumPz(0),SumE(0); 
  for(unsigned int i=0; i<tracks.size(); i++){
    SumPx += tracks.at(i)->px();              
    SumPy += tracks.at(i)->py();              
    SumPz += tracks.at(i)->pz();              
    SumE += sqrt(pow(tracks.at(i)->p(),2)+pow(PMH.Get_piMass(),2)); 
  }
  return sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2));
}

template <typename T>  bool ThreeProngInputSelector_Step1::cmpPt(const T & a, const T & b){
  return a->pt() > b->pt();
}


//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step1);




