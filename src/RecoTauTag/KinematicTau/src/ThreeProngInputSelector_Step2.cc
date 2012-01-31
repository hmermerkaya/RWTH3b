#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step2.h"

ThreeProngInputSelector_Step2::ThreeProngInputSelector_Step2(const edm::ParameterSet & iConfig):
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
    
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

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
std::vector<reco::TransientTrack> ThreeProngInputSelector_Step2::convToTransTrck(std::vector<reco::TrackRef> &input) {
	std::vector<reco::TransientTrack> transTrkVct;
	for (std::vector<reco::TrackRef>::iterator iter = input.begin(); iter != input.end(); ++iter) {
		transTrkVct.push_back( transTrackBuilder_->build( *iter ) );
	}
	return transTrkVct;
}
bool ThreeProngInputSelector_Step2::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx) {
	if (trkVct.size() < 2) {
		LogTrace("ThreeProngInputSelector_Step2") << "Can't check SecVertex: Only " << trkVct.size() << " Tracks.";
		return false;
	}else{
		bool useAdaptive = false;
		if (useAdaptive) {
			AdaptiveVertexFitter avf;//more robust?
			avf.setWeightThreshold(0.1);//weight per track. allow almost every fit, else --> exception
			try{
				transVtx = avf.vertex(trkVct);//AdaptiveVertexFitter
			}catch(...) {
                edm::LogWarning("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}else{
			KalmanVertexFitter kvf(true);
			try{
				transVtx = kvf.vertex(trkVct);//KalmanVertexFitter
			}catch(...) {
                edm::LogWarning("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}
		if (!transVtx.isValid()) LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::checkSecVtx: Secondary vertex not valid.";
		if (!useAdaptive) {
			if (!transVtx.hasRefittedTracks()) {
				LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::checkSecVtx: Secondary has 0 refitted tracks.";
				return false;
			}else if (transVtx.refittedTracks().size() != trkVct.size()) {
				LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::checkSecVtx: Secondary has only " << transVtx.refittedTracks().size() << " refitted of " << trkVct.size() << " initial tracks.";
				return false;
			}
		}
		
		return transVtx.isValid();
	}
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
            } else {
                ++pfiter;
            }
        }
        if (selected.size() >= minTau_) {
            LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: " << selected.size() << " tau candidate(s) reconstructed.";
            cntFound_++;
            return true;
        } else {
            LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Warning! Only " << selected.size() << " tau candidate(s) reconstructed. Skip Event.";
            return false;
        }
    }
    LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Unable to calculate a new primary vertex. No tau candidates reconstructed.";
    return false;
}
bool ThreeProngInputSelector_Step2::choose3bestTracks(std::vector<reco::TrackRefVector> & selected, std::vector<std::vector<reco::TrackRef> > combis, const reco::Vertex & pVtx) {
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
	if (combis.size()<1) {
		LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::choose3bestTracks: No combi survived.";
		return false;
	}
	
    reco::TrackRefVector tmpvec;
    
	if (combis.size()>1) {//chose the pair with smallest vertex rotation needed!!!
		sort(movements.begin(), movements.end(), pairSecond<int, double>);
		LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::choose3bestTracks:Too much combis (" << combis.size() << ") left. Take one with smallest vtx correction movement (" << movements.front().second << " [sigma]), second best is (" << movements.at(1).second << " [sigma]).";
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
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step2);
