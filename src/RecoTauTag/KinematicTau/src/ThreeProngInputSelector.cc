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
	produces<int>("threeProngFlag");//0=invalid, 1=valid
	produces<InputTrackCollection>("InputTracks");//save collection of vector<reco::CandidateRef> for each tau cand
    produces<InputTauCollection>("InputTauRefs");//needed to fill in unfit KinematicParticle later on
    produces<reco::VertexCollection>("primVtx");//has to be vector. save one of length one
}


ThreeProngInputSelector::~ThreeProngInputSelector(){
}



// ------------ method called on each new Event  ------------
bool ThreeProngInputSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	cnt_++;
	iEvent_ = &iEvent;
    
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

	std::auto_ptr<InputTrackCollection> selected_ = std::auto_ptr<InputTrackCollection >(new InputTrackCollection);
	InputTrackCollection & selected = * selected_;
    std::auto_ptr<InputTauCollection> PFTauRef_ = std::auto_ptr<InputTauCollection>(new InputTauCollection);
	InputTauCollection & PFTauRefs = * PFTauRef_;
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
	iEvent_->put(flagPtr,"threeProngFlag");
	
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
	printf("--> [ThreeProngInputSelector] found at least %i 3-prongs per event. Efficiency: %d/%d = %.2f%%\n", minTau_, cntFound_, cnt_, ratio*100.0);
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
		printf("Can't check SecVertex: Only %i Tracks.", trkVct.size());
		return false;
	}else{
		bool useAdaptive = false;
		if(useAdaptive){
			AdaptiveVertexFitter avf;//more robust?
			avf.setWeightThreshold(0.1);//weight per track. allow almost every fit, else --> exception
			try{
				transVtx = avf.vertex(trkVct);//AdaptiveVertexFitter
			}catch(...){
				printf("ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.\n");
				return false;
			}
		}else{
			KalmanVertexFitter kvf(true);
			try{
				transVtx = kvf.vertex(trkVct);//KalmanVertexFitter
			}catch(...){
				printf("ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.\n");
				return false;
			}
		}
		if(!transVtx.isValid()) printf("ThreeProngTauCreator::checkSecVtx: Secondary vertex not valid.\n");
		if(!useAdaptive){
			if(!transVtx.hasRefittedTracks()){
				LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::checkSecVtx: Secondary has 0 refitted tracks.";
				return false;
			}else if(transVtx.refittedTracks().size()!=trkVct.size()){
				LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::checkSecVtx: Secondary has only "<<transVtx.refittedTracks().size()<<" refitted of "<<trkVct.size()<<" initial tracks.";
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
			LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3Prongs: erased combi due to wrong charge sum. "<<combis.size()<<" combis left.";
			continue;
		}
		double massA1 = getInvariantMass(*iter, 0.140);
		if(massA1 > 2.0 || massA1 < 3*0.140){//soft upper value
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3Prongs: erased combi due to wrong mass. "<<combis.size()<<" combis left.";
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
        LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::createNewPrimVtx: No PrimaryVertexCollection found!";
        return false;
    }
    
    edm::InputTag beamSpotLabel = iConfig_.getParameter<edm::InputTag>("beamSpotLabel");
    reco::BeamSpot vertexBeamSpot;
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent_->getByLabel(beamSpotLabel,recoBeamSpotHandle);
    if (recoBeamSpotHandle.isValid()){
        vertexBeamSpot = *recoBeamSpotHandle;
    }else{
        edm::LogError("ThreeProngInputSelector") << "ThreeProngTauCreator::createNewPrimVtx: No beam spot available from EventSetup";
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
bool ThreeProngInputSelector::select(InputTrackCollection & selected, InputTauCollection & taurefs, reco::VertexCollection & primaryVertex){
    std::vector<std::vector<std::vector<reco::TrackRef> > > threeProngCombis;
    std::vector<reco::TrackRef> tautracks;
    
    edm::Handle<InputTrackCollection> inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
    
    if(!inputCollection.isValid()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::select: no InputTrackCollection found!";
        return false;
    }
    edm::Handle<InputTauCollection> inputTauCollection;
	iEvent_->getByLabel(selectedTauCandidatesTag_, inputTauCollection);
    taurefs = *inputTauCollection;
    
    if(!inputTauCollection.isValid()){
        LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::select: no InputTauCollection found!";
        return false;
    }
    
    
    for(InputTrackCollection::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
        threeProngCombis.push_back(choose3Prongs(input));
	}
    std::vector<std::vector<std::vector<reco::TrackRef> > >::const_iterator candidates;
    std::vector<std::vector<reco::TrackRef> > ::const_iterator triplets;
    std::vector<reco::TrackRef>::const_iterator tracks, trackrefs;
        
    for(candidates = threeProngCombis.begin(); candidates != threeProngCombis.end(); ++candidates){
        for(triplets = candidates->begin(); triplets != candidates->end(); ++triplets){
            for(tracks = triplets->begin(); tracks != triplets->end(); ++tracks){
                bool exists = false;
                for(trackrefs = tautracks.begin(); trackrefs != tautracks.end(); ++trackrefs){
                    if(*tracks == *trackrefs){
                        exists = true;
                        break;
                    }
                }
                if(!exists) tautracks.push_back(*tracks);
            }
        }
    }
        
    bool newvtxsuccess = createNewPrimVtx(primaryVertex, tautracks);
    if(newvtxsuccess) {
        std::vector<std::vector<std::vector<reco::TrackRef> > >::iterator candidates;
        InputTauCollection::iterator pfiter = taurefs.begin();
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
bool ThreeProngInputSelector::choose3bestTracks(InputTrackCollection & selected, std::vector<std::vector<reco::TrackRef> > combis, const reco::Vertex & pVtx){
	std::vector<std::pair<int,float> > chi2s;
	std::vector<std::pair<int,double> > movements;
	unsigned index=0;
	for (std::vector<std::vector<reco::TrackRef> >::iterator iter=combis.begin(); iter!=combis.end();) {
		TransientVertex tmpVtx;
		std::vector<reco::TransientTrack> trks = convToTransTrck(*iter);
		if(!checkSecVtx(trks, tmpVtx)){
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3bestTracks: Erased combi due to bad vertex. "<<combis.size()<<" combis left.";
			continue;
		}
        double massA1 = getInvariantMass(*iter, 0.140);
		TLorentzVector lorentzA1 = getSumTLorentzVec(*iter, massA1);
		VertexRotation vtxC(lorentzA1);
		double theta0;
		TVector3 tauFlghtDir;
		reco::Vertex pvTemp = pVtx;//do not modify original pv here
		vtxC.tryCorrection(pvTemp, tmpVtx, theta0, tauFlghtDir, false);//do not force rotation, only rotate within errors
		if(vtxC.isValid()) movements.push_back(std::make_pair(index,vtxC.movement()));
		else{
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3bestTracks: Erased combi due to bad vertex correction. "<<combis.size()<<" combis left.";
			continue;
		}
		
		chi2s.push_back(std::make_pair(index,tmpVtx.normalisedChiSquared()));
		
		++index;
		++iter;//only moved if nothing was deleted
	}
	if (combis.size()<1){
		LogTrace("ThreeProngTauCreator::choose3bestTracks: No combi survived.");
		return false;
	}
	
    reco::TrackRefVector tmpvec;
    
	if (combis.size()>1){//chose the pair with smallest vertex rotation needed!!!
		//		sort(chi2s.begin(), chi2s.end(), pairSecond<int, float>);
		//		LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3bestTracks:Too much combis ("<<combis.size()<<") left. Take best common vtx (chi2/ndf = "<<chi2s.at(0).second<<"), second best is (chi2/ndf = "<<chi2s.at(1).second<<").";
		//		unsigned int i = chi2s.front().first;
		//		input.assign( combis.at(i).begin(), combis.at(i).end() );
		
		sort(movements.begin(), movements.end(), pairSecond<int, double>);
		LogTrace("ThreeProngInputSelector")<<"ThreeProngTauCreator::choose3bestTracks:Too much combis ("<<combis.size()<<") left. Take one with smallest vtx correction movement ("<<movements.front().second<<" [sigma]), second best is ("<<movements.at(1).second<<" [sigma]).";
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




//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector);
