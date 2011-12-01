#include "RecoTauTag/KinematicTau/interface/InputTrackSelector.h"

InputTrackSelector::InputTrackSelector(const edm::ParameterSet& iConfig):
tauType_( iConfig.getUntrackedParameter<std::string>("tauType", "shrinkingCone") ),
minTracks_( iConfig.getParameter<unsigned int>("minTracks") ),//only tau candidates with more/equal than minTracks are selected
minTau_( iConfig.getUntrackedParameter<unsigned int>("minTau", 1) )//filter returns true if more/equal than minTau_ taus were selected
{
	produces<int>("flag");//0=invalid, 1=valid
	produces<std::vector<reco::TrackRefVector> >("InputTracks");//save collection of vector<reco::CandidateRef> for each tau cand
	produces<reco::PFTauRefVector>("InputTauRefs");//needed to fill in unfit KinematicParticle later on
}


InputTrackSelector::~InputTrackSelector(){
}



// ------------ method called on each new Event  ------------
bool InputTrackSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	cnt_++;
	iEvent_ = &iEvent;

	std::auto_ptr<std::vector<reco::TrackRefVector> > selected_ = std::auto_ptr<std::vector<reco::TrackRefVector> >(new std::vector<reco::TrackRefVector>);
	std::vector<reco::TrackRefVector> & selected = * selected_;
	std::auto_ptr<reco::PFTauRefVector> PFTauRef_ = std::auto_ptr<reco::PFTauRefVector>(new reco::PFTauRefVector);
	reco::PFTauRefVector & PFTauRef = * PFTauRef_;

	bool filterValue = select(selected, PFTauRef);

	iEvent_->put(PFTauRef_,"InputTauRefs");
	iEvent_->put(selected_,"InputTracks");

	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(filterValue) flag = 1;
	else flag = 0;
	iEvent_->put(flagPtr,"flag");
	
	return filterValue;
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
    edm::LogVerbatim("InputTrackSelector")<<"--> [InputTrackSelector] found at least "<<minTau_<<" tau candidate(s) of type "<<tauType_<<" (with at least "<<minTracks_<<" tracks) per event. Efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}
bool InputTrackSelector::select(std::vector<reco::TrackRefVector> & selected, reco::PFTauRefVector & PFTauRef){
	bool found = false;
	
	edm::Handle<reco::PFTauCollection> inputCollection;
	iEvent_->getByLabel(tauType_+"PFTauProducer", inputCollection);
	LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: no. of PFTaus in event = "<<inputCollection->size();
	for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
		reco::PFTauRef thePFTau(inputCollection, iPFTau);
		LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: test tau "<<iPFTau;
		reco::TrackRefVector tauDaughters = getPFTauDaughters(thePFTau);
		if(tauDaughters.size()>=minTracks_){
			selected.push_back(tauDaughters);
			PFTauRef.push_back(thePFTau);
		}else LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: only "<<tauDaughters.size()<<" tau daughter(s) found. Skip tau candidate.";
	}
	
	if(selected.size() >= minTau_){
		cntFound_++;
		found = true;
		LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: "<<selected.size()<<" tau candidate(s) reconstructed.";
	}else LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select:Warning: only "<<selected.size()<<" tau candidate(s) reconstructed. Skip Evt.";

	return found;
}
reco::TrackRefVector InputTrackSelector::getPFTauDaughters(reco::PFTauRef &PFTau){
	reco::TrackRefVector trkVct;
	const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands(); //candidates in signal cone 
	//isolationPFChargedHadrCands stores tracks in isolation/veto cone
	for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
		LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::getPFTauDaughters: PFTau daughter pt "<<iter->get()->pt()<<", eta "<<iter->get()->eta()<<", vtx("<<iter->get()->vx()<<","<<iter->get()->vy()<<","<<iter->get()->vz()<<")";
		if(iter->get()->trackRef().isNonnull()) trkVct.push_back( (*iter)->trackRef() );
	}
	
	return trkVct;
}


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);
