#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector.h"

ThreeProngInputSelector::ThreeProngInputSelector(const edm::ParameterSet& iConfig):
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauCandidates" ) ),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
selectedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "selectedTauCandidates" ) ),
minTau_( iConfig.getUntrackedParameter<unsigned int>("minTau", 1) ),//filter returns true if more/equal than minTau_ taus were selected
minVtxTracks_(iConfig.getUntrackedParameter("minVtxTracks", int(3))),
maxChi2ndf_(iConfig.getUntrackedParameter("maxChi2ndf", double(10.0)))
{
    iConfig_ = iConfig;
	produces<int>("flag");//0=invalid, 1=valid
	produces<std::vector<reco::TrackRefVector> >("InputTracks");//save collection of vector<reco::CandidateRef> for each tau cand
    produces<reco::PFTauRefVector>("InputTauRefs");//needed to fill in unfit KinematicParticle later on
    produces<reco::VertexCollection>("primVtx");//has to be vector. save one of length one
}


ThreeProngInputSelector::~ThreeProngInputSelector(){
}



// ------------ method called on each new Event  ------------
bool ThreeProngInputSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	cnt_++;
	iEvent_ = &iEvent;
    
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

	std::auto_ptr<std::vector<reco::TrackRefVector> > selected_ = std::auto_ptr<std::vector<reco::TrackRefVector> >(new std::vector<reco::TrackRefVector>);
	std::vector<reco::TrackRefVector> & selected = * selected_;
    std::auto_ptr<reco::PFTauRefVector> PFTauRef_ = std::auto_ptr<reco::PFTauRefVector>(new reco::PFTauRefVector);
	reco::PFTauRefVector & PFTauRefs = * PFTauRef_;
    std::auto_ptr<reco::VertexCollection> primaryVertex_ = std::auto_ptr<reco::VertexCollection>(new reco::VertexCollection);
	reco::VertexCollection & primaryVertex = * primaryVertex_;
        
	bool filterValue = select(selected, PFTauRefs, primaryVertex);
    
	iEvent_->put(selected_,"InputTracks");
    iEvent_->put(PFTauRef_,"InputTauRefs");
	iEvent_->put(primaryVertex_,"primVtx");

	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(filterValue) flag = 1;
	else flag = 0;
	iEvent_->put(flagPtr,"flag");
	
	return filterValue;
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeProngInputSelector::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeProngInputSelector::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
    edm::LogVerbatim("ThreeProngInputSelector")<<"--> [ThreeProngInputSelector] found at least "<<minTau_<<" 3-prongs per event. Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}
bool ThreeProngInputSelector::sumCharge(std::vector<reco::TrackRef> &input){
	int sum = abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge());
	if(sum == 1) return true;
	return false;
}
std::vector<reco::TransientTrack> ThreeProngInputSelector::convToTransTrck(std::vector<reco::TrackRef> &input){
	std::vector<reco::TransientTrack> transTrkVct;
	for (std::vector<reco::TrackRef>::iterator iter=input.begin(); iter!=input.end(); ++iter) {
		transTrkVct.push_back( transTrackBuilder_->build( *iter ) );
	}
	return transTrkVct;
}
bool ThreeProngInputSelector::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx){
	if(trkVct.size()<2){
		LogTrace("ThreeProngInputSelector")<<"Can't check SecVertex: Only "<<trkVct.size()<<" Tracks.";
		return false;
	}else{
		bool useAdaptive = false;
		if(useAdaptive){
			AdaptiveVertexFitter avf;//more robust?
			avf.setWeightThreshold(0.1);//weight per track. allow almost every fit, else --> exception
			try{
				transVtx = avf.vertex(trkVct);//AdaptiveVertexFitter
			}catch(...){
                edm::LogWarning("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}else{
			KalmanVertexFitter kvf(true);
			try{
				transVtx = kvf.vertex(trkVct);//KalmanVertexFitter
			}catch(...){
                edm::LogWarning("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}
		if(!transVtx.isValid()) LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkSecVtx: Secondary vertex not valid.";
		if(!useAdaptive){
			if(!transVtx.hasRefittedTracks()){
				LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkSecVtx: Secondary has 0 refitted tracks.";
				return false;
			}else if(transVtx.refittedTracks().size()!=trkVct.size()){
				LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkSecVtx: Secondary has only "<<transVtx.refittedTracks().size()<<" refitted of "<<trkVct.size()<<" initial tracks.";
				return false;
			}
		}
		
		return transVtx.isValid();
	}
}
template <typename T> std::vector<std::vector<T> > ThreeProngInputSelector::permuteCombinations(const std::vector<T> &vect){
	std::vector<std::vector<T> > combis;
	typename std::vector<T>::const_iterator iter1, iter2, iter3;
	for (iter1 = vect.begin(); iter1!=vect.end(); ++iter1) {
		iter2 = iter1;
		++iter2;
		for (; iter2!=vect.end(); ++iter2) {
			iter3 = iter2;
			++iter3;
			for (; iter3!=vect.end(); ++iter3) {
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
std::vector<std::vector<reco::TrackRef> > ThreeProngInputSelector::choose3Prongs(std::vector<reco::TrackRef> &input){
    sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);
	std::vector<std::vector<reco::TrackRef> > combis = permuteCombinations(input);

	for (std::vector<std::vector<reco::TrackRef> >::iterator iter=combis.begin(); iter!=combis.end();) {
		if(!sumCharge(*iter)){
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::choose3Prongs: erased combi due to wrong charge sum. "<<combis.size()<<" combis left.";
			continue;
		}
		double massA1 = getInvariantMass(*iter, 0.140);
		if(massA1 > 2.0 || massA1 < 3*0.140){//soft upper value
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::choose3Prongs: erased combi due to wrong mass. "<<combis.size()<<" combis left.";
			continue;
		}
        ++iter;
    }
    
    return combis;
}
bool ThreeProngInputSelector::createNewPrimVtx(reco::VertexCollection & primaryVertex, const std::vector<reco::TrackRef> & tautracks){
    edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);

    if(!primVtxs.isValid()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::createNewPrimVtx: No PrimaryVertexCollection found!";
        return false;
    }
    
    edm::InputTag beamSpotLabel = iConfig_.getParameter<edm::InputTag>("beamSpotLabel");
    reco::BeamSpot vertexBeamSpot;
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent_->getByLabel(beamSpotLabel,recoBeamSpotHandle);
    if (recoBeamSpotHandle.isValid()){
        vertexBeamSpot = *recoBeamSpotHandle;
    }else{
        edm::LogError("ThreeProngInputSelector") << "ThreeProngInputSelector::createNewPrimVtx: No beam spot available from EventSetup";
    }
    
    PrimaryVertexProducerAlgorithm PVPA(iConfig_);
    std::vector<reco::TrackRef> vtxtracks;
    for (reco::VertexCollection::const_iterator iter = primVtxs->begin(); iter != primVtxs->end(); ++iter) {
        for (reco::Vertex::trackRef_iterator tracks = iter->tracks_begin(); tracks != iter->tracks_end(); ++tracks) {
            bool exclude = false;
            for(std::vector<reco::TrackRef>::const_iterator tautrackrefs = tautracks.begin(); tautrackrefs != tautracks.end(); ++tautrackrefs){
                if((*tracks).castTo<reco::TrackRef>() == *tautrackrefs){
                    exclude = true;
                    break;
                }
            }
            if(!exclude) vtxtracks.push_back((*tracks).castTo<reco::TrackRef>());
        }
    }
    std::vector<TransientVertex> newvertices = PVPA.vertices(convToTransTrck(vtxtracks), vertexBeamSpot);
    return checkPrimVtx(primaryVertex, newvertices);
}
bool ThreeProngInputSelector::select(std::vector<reco::TrackRefVector> & selected, reco::PFTauRefVector & taurefs, reco::VertexCollection & primaryVertex){
    std::vector<std::vector<std::vector<reco::TrackRef> > > threeProngCombis;
    std::vector<reco::TrackRef> tautracks;
    
	//load input collection from InputTrackSelector
	//each tau candidate stores its signal cone tracks
    edm::Handle<std::vector<reco::TrackRefVector> > inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
    
    if(!inputCollection.isValid()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: no std::vector<reco::TrackRefVector> found!";
        return false;
    }
    edm::Handle<reco::PFTauRefVector> inputTauCollection;
	iEvent_->getByLabel(selectedTauCandidatesTag_, inputTauCollection);
    taurefs = *inputTauCollection;
    
    if(!inputTauCollection.isValid()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: no inputTauCollection found!";
        return false;
    }
	if(inputCollection->size() != inputTauCollection->size()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<selectedTauCandidatesTag_.label()<<"("<<inputTauCollection->size()<<")";
		return false;
	}
	
    
    
    for(std::vector<reco::TrackRefVector>::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
        threeProngCombis.push_back(choose3Prongs(input));//make combinations of exact 3 tracks, choose3Prongs can return an empty vector which will be deleted later on (together with its tauRef)
	}

	//save all signal cone tracks of all candidates ONCE into tautracks without duplicates
	//keep duplicates in duplicateTracks and delete combi if all three tracks are also in ONE other combi (of another candidate)
	std::vector<std::vector<std::vector<reco::TrackRef> > >::iterator candidates, passedCandidates;
    std::vector<std::vector<reco::TrackRef> > ::iterator triplets, passedTriplets;
    std::vector<reco::TrackRef>::const_iterator tracks, trackrefs;
    for(candidates = threeProngCombis.begin(); candidates != threeProngCombis.end(); ++candidates){
        for(triplets = candidates->begin(); triplets != candidates->end();){
			bool erasedTriplet = false;
			std::vector<reco::TrackRef> duplicateTracks;
            for(tracks = triplets->begin(); tracks != triplets->end(); ++tracks){
                bool exists = false;
                for(trackrefs = tautracks.begin(); trackrefs != tautracks.end(); ++trackrefs){
                    if(*tracks == *trackrefs){
                        exists = true;
                        break;
                    }
                }
                if(!exists) tautracks.push_back(*tracks);
				else duplicateTracks.push_back(*tracks);
            }
			if(duplicateTracks.size()==3) erasedTriplet = removeDuplicateTriplets(duplicateTracks, threeProngCombis, candidates, triplets);
			if(!erasedTriplet) ++triplets;
        }
    }
        
    bool newvtxsuccess = createNewPrimVtx(primaryVertex, tautracks);
    if(newvtxsuccess) {
        reco::PFTauRefVector::iterator pfiter = taurefs.begin();
        for (candidates = threeProngCombis.begin(); candidates != threeProngCombis.end(); ++candidates) {
            bool success = choose3bestTracks(selected, *candidates, primaryVertex.front());
            if (!success) {
                pfiter = taurefs.erase(pfiter);
            } else {
                ++pfiter;
            }
        }
        if (selected.size() >= minTau_) {
            LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: "<<selected.size()<<" tau candidate(s) reconstructed.";
            cntFound_++;
            return true;
        } else {
            LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: Warning! Only "<<selected.size()<<" tau candidate(s) reconstructed. Skip Event.";
            return false;
        }
    }
    LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::select: Unable to calculate a new primary vertex. No tau candidates reconstructed.";
    return false;
}
bool ThreeProngInputSelector::checkPrimVtx(reco::VertexCollection & primaryVertex, const std::vector<TransientVertex> & newvertices){
    std::vector<reco::Vertex> vtx;
    for (vector<TransientVertex>::const_iterator iv = newvertices.begin(); iv != newvertices.end(); ++iv) {
        reco::Vertex v = *iv;
        if(v.tracksSize() >= minVtxTracks_ && v.normalizedChi2() <= maxChi2ndf_) vtx.push_back(v);
    }
	if(vtx.size()<1){
		LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::checkPrimVtx: No valid primary vertex found. Skip event "<<iEvent_->id().event()<<".";
		return false;
	}
	if(vtx.size()>1){
		sort(vtx.begin(), vtx.end(), cmpNormalizedChi2<reco::Vertex>);
	}
	primaryVertex.push_back(vtx.front());
	return true;
}
bool ThreeProngInputSelector::choose3bestTracks(std::vector<reco::TrackRefVector> & selected, std::vector<std::vector<reco::TrackRef> > combis, const reco::Vertex & pVtx){
	std::vector<std::pair<int,double> > movements;
	unsigned index=0;
	for (std::vector<std::vector<reco::TrackRef> >::iterator iter=combis.begin(); iter!=combis.end();) {
		TransientVertex tmpVtx;
		std::vector<reco::TransientTrack> trks = convToTransTrck(*iter);
		if(!checkSecVtx(trks, tmpVtx)){
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::choose3bestTracks: Erased combi due to bad vertex. "<<combis.size()<<" combis left.";
			continue;
		}
        double massA1 = getInvariantMass(*iter, 0.140);
		TLorentzVector lorentzA1 = getSumTLorentzVec(*iter, massA1);
		VertexRotation vtxC(lorentzA1);
		double theta0;
		TVector3 tauFlghtDir;
		reco::Vertex pvTemp = pVtx;//do not modify original pv here
		double significance = vtxC.rotatePV(pvTemp, tmpVtx, theta0, tauFlghtDir);
		movements.push_back(std::make_pair(index, significance));
		
		++index;
		++iter;//only moved if nothing was deleted
	}
	if (combis.size()<1){
		LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::choose3bestTracks: No combi survived.";
		return false;
	}
	
    reco::TrackRefVector tmpvec;
    
	if (combis.size()>1){//chose the pair with smallest vertex rotation needed!!!
		sort(movements.begin(), movements.end(), pairSecond<int, double>);
		LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::choose3bestTracks:Too much combis ("<<combis.size()<<") left. Take one with smallest vtx correction movement ("<<movements.front().second<<" [sigma]), second best is ("<<movements.at(1).second<<" [sigma]).";
		unsigned int i = movements.front().first;
        for (std::vector<reco::TrackRef>::const_iterator iter = combis.at(i).begin(); iter != combis.at(i).end(); ++iter) {
            tmpvec.push_back(*iter);
        }
		
	} else {
        for (std::vector<reco::TrackRef>::const_iterator iter = combis.front().begin(); iter != combis.front().end(); ++iter) {
            tmpvec.push_back(*iter);
        }
    }
    
    selected.push_back(tmpvec);
    
	return true;
}
bool ThreeProngInputSelector::removeDuplicateTriplets(const std::vector<reco::TrackRef> & duplicateTracks, 
													  std::vector<std::vector<std::vector<reco::TrackRef> > > & threeProngCombis, 
													  std::vector<std::vector<std::vector<reco::TrackRef> > >::iterator & candidates, 
													  std::vector<std::vector<reco::TrackRef> > ::iterator & triplets){
	//check on all already tested combis if all three tracks of this duplicate triplet belong to one COMMON other triplet
	//so that the whole triplet is equal
	//do not delete triplet if tracks belong to different other combis!

	std::vector<std::vector<std::vector<reco::TrackRef> > >::const_iterator passedCandidates;
    std::vector<std::vector<reco::TrackRef> > ::const_iterator passedTriplets;
    std::vector<reco::TrackRef>::const_iterator tracks, duplicate;	
	for(passedCandidates = threeProngCombis.begin(); passedCandidates != threeProngCombis.end(); ++passedCandidates){
		for(passedTriplets = passedCandidates->begin(); passedTriplets != passedCandidates->end(); ++passedTriplets){
			if(passedTriplets==triplets){//test only already passed triplets
				return false;
			}
			unsigned int cntDuplicate = 0;
			for(tracks = passedTriplets->begin(); tracks != passedTriplets->end(); ++tracks){
				for(duplicate = duplicateTracks.begin(); duplicate != duplicateTracks.end(); ++duplicate){
					if(*tracks == *duplicate){
						cntDuplicate++;
					}
				}
			}
			if(cntDuplicate==3){//All 3 tracks are repeated in one COMMON other triplet. Therefore delete this triplet.
				LogTrace("ThreeProngInputSelector")<<"ThreeProngInputSelector::removeDuplicateTriplets: Delete duplicate triplet!";
				triplets = candidates->erase(triplets);
				return true;
			}
		}
		if(passedCandidates==candidates){//test only already passed candidates
			return false;
		}
	}
    edm::LogError("ThreeProngInputSelector")<<"ThreeProngInputSelector::removeDuplicateTriplets: One should never see this!";
	return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector);
