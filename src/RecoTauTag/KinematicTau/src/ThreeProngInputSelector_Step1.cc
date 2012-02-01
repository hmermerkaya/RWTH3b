#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step1.h"

ThreeProngInputSelector_Step1::ThreeProngInputSelector_Step1(const edm::ParameterSet& iConfig):
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
	
    bool filterValue = select(nonTauTracks, threeProngCombis);
    
	iEvent_->put(pNonTauTracks, "NonTauTracks");
	iEvent_->put(pThreeProngCombis, "ThreeProngCombinations");

	std::auto_ptr<int> pFlag = std::auto_ptr<int>(new int);
	int & flag = * pFlag;
	if (filterValue) {
        flag = 1;
    } else {
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
bool ThreeProngInputSelector_Step1::sumCharge(const std::vector<reco::TrackRef> & input) {
	int sum = abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge());
	if (sum == 1) return true;
	return false;
}
template <typename T> std::vector<std::vector<T> > ThreeProngInputSelector_Step1::permuteCombinations(const std::vector<T> & vect) {
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
vVTrackRef ThreeProngInputSelector_Step1::choose3Prongs(std::vector<reco::TrackRef> & input) {
    sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);
	vVTrackRef combis = permuteCombinations(input);

	for (vVTrackRef::iterator iter = combis.begin(); iter != combis.end();) {
		if (!sumCharge(*iter)) {
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong charge sum. " << combis.size() << " combis left.";
			continue;
		}
		double massA1 = getInvariantMass(*iter, 0.140);
		if (massA1 > 2.0 || massA1 < 3*0.140) { //soft upper value
			iter = combis.erase(iter);
			LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::choose3Prongs: erased combi due to wrong mass. " << combis.size() << " combis left.";
			continue;
		}
        ++iter;
    }
    
    return combis;
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
                } else {
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
    } else {
        LogTrace("ThreeProngInputSelector_Step1") << "ThreeProngInputSelector_Step1::select: Warning! Only " << threeProngCombis.size() << " tau candidate(s) reconstructed. Skip Event.";
        return false;
    }
    return false;
}
bool ThreeProngInputSelector_Step1::removeDuplicateTriplets(const std::vector<reco::TrackRef> & duplicateTracks, vVVTrackRef & threeProngCombis, vVVTrackRef::iterator & candidates, vVTrackRef::iterator & triplets) {
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

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step1);
