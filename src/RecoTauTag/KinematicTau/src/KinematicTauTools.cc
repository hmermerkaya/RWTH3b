#include  "RecoTauTag/KinematicTau/interface/KinematicTauTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


const double KinematicTauTools::piMass = 0.13957018;
const double KinematicTauTools::tauMass = 1.77682;

KinematicTauTools::KinematicTauTools(){
}

KinematicTauTools::~KinematicTauTools(){
}

/////////////////////////////////////////////////////////////////////////////
// Compute the charge of the triplet 
bool KinematicTauTools::sumCharge(const std::vector<reco::TrackRef> & input){
  return (1== abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge()));
}


/////////////////////////////////////////////////////////////////////////////
// run all permutations
template <typename T> std::vector<std::vector<T> > KinematicTauTools::permuteCombinations(const std::vector<T> & vect){
  std::vector<std::vector<T> > combis;
  typename std::vector<T>::const_iterator iter1, iter2, iter3;
  for (iter1 = vect.begin(); iter1 != vect.end(); ++iter1) {
    iter2 = iter1;
    ++iter2;
    for (; iter2 != vect.end(); ++iter2) {
      iter3 = iter2;
      ++iter3;
      for (; iter3 != vect.end(); ++iter3) {
	std::vector<T> newCombi;
	newCombi.push_back(*iter1);
	newCombi.push_back(*iter2);
	newCombi.push_back(*iter3);
	combis.push_back(newCombi);
      }
    }
  }
  return combis;
}


/////////////////////////////////////////////////////////////////////////////
// Select 3 prong based on mass and charge
vVTrackRef KinematicTauTools::choose3Prongs(std::vector<reco::TrackRef> & input){
  sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);
  vVTrackRef combis = permuteCombinations(input);
  for (vVTrackRef::iterator iter = combis.begin(); iter != combis.end();) {
    if (!sumCharge(*iter)) {
      iter = combis.erase(iter);
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong charge sum. " << combis.size() << " combis left.";
      continue;
    }
    double massA1 = getInvariantMass(*iter);
    if (massA1 > 2.0 || massA1 < 3*Get_piMass()) { //soft upper value                                                     
      iter = combis.erase(iter);
      LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong mass. " << combis.size() << " combis left.";
      continue;
    }
    ++iter;
  }
  return combis;
}

/////////////////////////////////////////////////////////////////////////////
// remove dubplicate triplets
bool KinematicTauTools::removeDuplicateTriplets(const std::vector<reco::TrackRef> & duplicateTracks, vVVTrackRef & threeProngCombis, vVVTrackRef::iterator & candidates, vVTrackRef::iterator & triplets){
  //check on all already tested combis if all three tracks of this duplicate triplet belong to one COMMON other triplet                                                  
  //so that the whole triplet is equal 
  //do not delete triplet if tracks belong to different other combis!  
  std::vector<reco::TrackRef>::const_iterator tracks, duplicate;
  for (vVVTrackRef::const_iterator passedCandidates = threeProngCombis.begin(); passedCandidates != threeProngCombis.end(); ++passedCandidates) {
    for (vVTrackRef::const_iterator passedTriplets = passedCandidates->begin(); passedTriplets != passedCandidates->end(); ++passedTriplets) {
      if (passedTriplets == triplets) { //test only already passed triplets                                                                                                                        
	return false;
      }
      unsigned int cntDuplicate = 0;
      for (tracks = passedTriplets->begin(); tracks != passedTriplets->end(); ++tracks) {
	for (duplicate = duplicateTracks.begin(); duplicate != duplicateTracks.end(); ++duplicate) {
	  if (*tracks == *duplicate) {
	    cntDuplicate++;
	  }
	}
      }
      if (cntDuplicate == 3) { //All 3 tracks are repeated in one COMMON other triplet. Therefore delete this triplet.                                                                             
	LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::removeDuplicateTriplets: Delete duplicate triplet!";
	triplets = candidates->erase(triplets);
	return true;
      }
    }
    if (passedCandidates == candidates) { //test only already passed candidates                                                                                                                          
      return false;
    }
  }
  edm::LogError("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::removeDuplicateTriplets: One should never see this!";
  return false;
}

/////////////////////////////////////////////////////////////////////////////
// Check for secondary vertex
bool KinematicTauTools::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx){
  if(trkVct.size()<2){
    LogTrace("ThreeProngTauCreator")<<"Can't check SecVertex: Only "<<trkVct.size()<<" Tracks.";
    return false;
  }else{
    bool useAdaptive = false;
    if(useAdaptive){
      AdaptiveVertexFitter avf;
      avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception                                                                                                                          
      try{
	transVtx = avf.vertex(trkVct); //AdaptiveVertexFitter                                                                                                                                                        
      }catch(...){
	LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
	return false;
      }
    }else{
      KalmanVertexFitter kvf(true);
      try{
	transVtx = kvf.vertex(trkVct); //KalmanVertexFitter                                                                                                                                                          
      }catch(...){
	LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
	return false;
      }
    }
    if(!transVtx.isValid()) LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex not valid.";
    if(!useAdaptive){
      if(!transVtx.hasRefittedTracks()){
	LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has 0 refitted tracks.";
	return false;
      }else if(transVtx.refittedTracks().size()!=trkVct.size()){
	LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has only "<<transVtx.refittedTracks().size()<<" refitted of "<<trkVct.size()<<" initial tracks.";
	return false;
      }
    }

    return transVtx.isValid();
  }
}


std::vector<reco::TransientTrack> KinematicTauTools::convToTransTrck(std::vector<reco::TrackRef> &input){
  //TransientTrackBuilder delivers reco::TransientTrack
  std::vector<reco::TransientTrack> transTrkVct;
  for (std::vector<reco::TrackRef>::const_iterator iter=input.begin(); iter!=input.end(); ++iter) {
    transTrkVct.push_back( transientTrackBuilder_->build( *iter ) );
  }
  return transTrkVct;
}


/////////////////////////////////////////////////////////////////////////////
// compute the invariant mass
template <typename T> double KinematicTauTools::getInvariantMass(const T & tracks){//if second argument empty default pion is supposed
  double SumPx(0),SumPy(0),SumPz(0),SumE(0);
  for(unsigned int i=0; i<tracks.size(); i++){
    SumPx += tracks.at(i)->px();
    SumPy += tracks.at(i)->py();
    SumPz += tracks.at(i)->pz();
    SumE += sqrt(pow(tracks.at(i)->p(),2)+pow(Get_piMass(),2));
  }
  return sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2));
}

/////////////////////////////////////////////////////////////////////////////
// comparision for pt
template <typename T>  bool KinematicTauTools::cmpPt(const T & a, const T & b){
  return a->pt() > b->pt();
}


double KinematicTauTools::VertexRotationAndSignificance(const std::vector<reco::TrackRef> &input,TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks,reco::Vertex &pVtx,TLorentzVector &lorentzA1, TVector3 &tauFlghtDir, double &theta0, double &thetaMax){
  double massA1 = getInvariantMass(input);
  lorentzA1 = getSumTLorentzVec(input, massA1);
  VertexRotation vtxC(lorentzA1);
  thetaMax=fabs(vtxC.calcThetaMax());
  return vtxC.rotatePV(pVtx,tmpVtx,theta0, tauFlghtDir);
}


bool KinematicTauTools::choose3bestTracks(std::vector<reco::TrackRef> & input, reco::Vertex & pVtx){
  std::vector<std::vector<reco::TrackRef> > combis=KinematicTauTools::choose3Prongs(input);
  std::vector<reco::TrackRefVector> selected;
  return KinematicTauTools::choose3bestTracks(selected,combis,pVtx);
}


bool KinematicTauTools::choose3bestTracks(std::vector<reco::TrackRefVector> & selected, std::vector<std::vector<reco::TrackRef> > combis, const reco::Vertex &pVtx) {
  std::vector<std::pair<int,double> > movements;
  unsigned index = 0;
  for (std::vector<std::vector<reco::TrackRef> >::iterator iter = combis.begin(); iter != combis.end();) {
    TransientVertex tmpVtx;
    std::vector<reco::TransientTrack> trks = convToTransTrck(*iter);
    if (!checkSecVtx(trks, tmpVtx)) {
      iter = combis.erase(iter);
      LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::choose3bestTracks: Erased combi due to bad vertex. " << combis.size() << " combis left.";
      continue;
    }
    TLorentzVector lorentzA1;
    double theta0,thetaMax;
    TVector3 tauFlghtDir;
    reco::Vertex pvTemp = pVtx;//do not modify original pv here
    double significance = VertexRotationAndSignificance(*iter,tmpVtx,trks,pvTemp,lorentzA1,tauFlghtDir,theta0,thetaMax);

    edm::LogInfo("ThreeProngInputSelector_Step2")<<"ThreeProngInputSelector_Step2::choose3bestTracks Original method significance " << significance
						 << " PVertexFit and Rotate (" 
						 <<  pVtx.position().x()
						 << "," <<  pVtx.position().y()
						 << "," <<  pVtx.position().z() << ")"
						 << " SV (" <<  tmpVtx.position().x()
						 << "," <<  tmpVtx.position().y()
						 << "," <<  tmpVtx.position().z() << ")";

    movements.push_back(std::make_pair(index, significance));

    ++index;
    ++iter;//only moved if nothing was deleted                                                                                                                                                                                               
  }
  if (combis.size()<1) {
    LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::choose3bestTracks: No combi survived.";
    return false;
  }
  reco::TrackRefVector tmpvec;
  if (combis.size()>1) {//chose the pair with smallest vertex rotation needed!!!                                                                                                                                                             
    sort(movements.begin(), movements.end(), pairSecond<int, double>);
    LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::choose3bestTracks:Too much combis (" << combis.size() << ") left. Take one with smallest vtx correction movement ("
                                              << movements.front().second << " [sigma]), second best is (" << movements.at(1).second << " [sigma]).";
    unsigned int i = movements.front().first;
    for (std::vector<reco::TrackRef>::const_iterator iter = combis.at(i).begin(); iter != combis.at(i).end(); ++iter) {
      tmpvec.push_back(*iter);
    }
  }
  else {
    for (std::vector<reco::TrackRef>::const_iterator iter = combis.front().begin(); iter != combis.front().end(); ++iter) {
      tmpvec.push_back(*iter);
    }
  }
  selected.push_back(tmpvec);
  return true;
}





template <class T> TLorentzVector KinematicTauTools::getSumTLorentzVec(const T& tracks, const double massConstraint){
  double sumPx=0, sumPy=0, sumPz=0;
  for(unsigned int i=0; i!= tracks.size(); i++){
    sumPx += tracks.at(i)->px();
    sumPy += tracks.at(i)->py();
    sumPz += tracks.at(i)->pz();
  }
  TLorentzVector lorentz;
  lorentz.SetXYZM(sumPx, sumPy, sumPz, massConstraint);
  return lorentz;
}

template <typename S, typename T> bool KinematicTauTools::pairSecond(const std::pair<S, T> &a, const std::pair<S, T> &b){
  return a.second < b.second;
}

template <typename T> bool KinematicTauTools::cmpChi2(const T &a, const T &b){
  return a->chi2() > b->chi2();
}

bool  KinematicTauTools::GetNonTauTracks(edm::Event *iEvent_,edm::InputTag &trackCollectionTag_,reco::TrackCollection &nonTauTracks, std::vector<reco::TrackRef> &tautracks){
  //load general track collection and substract tau tracks from it 
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

