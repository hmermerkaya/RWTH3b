#include "RecoTauTag/KinematicTau/interface/InputTrackSelector.h"

InputTrackSelector::InputTrackSelector(const edm::ParameterSet& iConfig):
tauType_( iConfig.getUntrackedParameter<std::string>("tauType", "fixedCone") ),
minTracks_( iConfig.getParameter<unsigned int>("minTracks") ),//only tau candidates with more/equal than minTracks are selected
minTau_( iConfig.getUntrackedParameter<unsigned int>("minTau", 1) ),//filter returns true if more/equal than minTau_ taus were selected
filterTaus_( iConfig.getUntrackedParameter<bool>("filterTaus", true) )//decide whether to pre-filter the pftaus, default is true
{
	produces<int>("flag");//0=invalid, 1=valid
	produces<InputTrackCollection>("InputTracks");//save collection of vector<reco::CandidateRef> for each tau cand
	produces<InputTauCollection>("InputTauRefs");//needed to fill in unfit KinematicParticle later on
}


InputTrackSelector::~InputTrackSelector(){
}



// ------------ method called on each new Event  ------------
bool InputTrackSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	cnt_++;
	iEvent_ = &iEvent;

	std::auto_ptr<InputTrackCollection> selected_ = std::auto_ptr<InputTrackCollection >(new InputTrackCollection);
	InputTrackCollection & selected = * selected_;
	std::auto_ptr<InputTauCollection> PFTauRef_ = std::auto_ptr<InputTauCollection>(new InputTauCollection);
	InputTauCollection & PFTauRef = * PFTauRef_;

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
	printf("--> [InputTrackSelector] found at least %i tau candidate (with at least %i tracks) per event. Efficiency: %d/%d = %.2f%%\n", minTau_, minTracks_, cntFound_, cnt_, ratio*100.0);
}
bool InputTrackSelector::select(InputTrackCollection & selected, InputTauCollection & PFTauRef){
	bool found = false;
	
	edm::Handle<reco::PFTauCollection> inputCollection;
	iEvent_->getByLabel(tauType_+"PFTauProducer", inputCollection);
	LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::select: no. of PFTaus in event = "<<inputCollection->size();
	for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
		reco::PFTauRef thePFTau(inputCollection, iPFTau);
		if(filterTaus_){
			if(!filterInput(thePFTau)) continue;//move into external module?
		}

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
bool InputTrackSelector::filterInput(reco::PFTauRef &tau){//use seperate filter module?
	bool filter = false;
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorByIsolation;
	iEvent_->getByLabel(tauType_+"PFTauDiscriminationByTrackIsolation",thePFTauDiscriminatorByIsolation); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorByLeadingTrackPtCut; 
	iEvent_->getByLabel(tauType_+"PFTauDiscriminationByLeadingTrackPtCut",thePFTauDiscriminatorByLeadingTrackPtCut); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorAgainstElectrons; 
	iEvent_->getByLabel(tauType_+"PFTauDiscriminationAgainstElectron",thePFTauDiscriminatorAgainstElectrons); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorAgainstMuons; 
	iEvent_->getByLabel(tauType_+"PFTauDiscriminationAgainstMuon",thePFTauDiscriminatorAgainstMuons);
	
//edmDumpEventContent /disk1/perchalla/data/CMSSW_3_1_2/KinTau/tau3piFromVBFH/AODSIMHLT_tau3piFromVBFH_145GeV.root | grep PFTau	
//	fixedConeHighEffPFTauDiscriminationAgainstElectron
//	fixedConeHighEffPFTauDiscriminationAgainstMuon
//	fixedConeHighEffPFTauDiscriminationByECALIsolation
//	fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion
//	fixedConeHighEffPFTauDiscriminationByIsolation
//	fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion
//	fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut
//	fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding
//	fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut //MinPtLeadingTrack = cms.double(5.0)
//	fixedConeHighEffPFTauDiscriminationByTrackIsolation
//	fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion
	
    if (
        ((*thePFTauDiscriminatorAgainstMuons)[tau] == 1 || (*thePFTauDiscriminatorAgainstElectrons)[tau] == 1)
        && ((*thePFTauDiscriminatorByIsolation)[tau] == 1 || (*thePFTauDiscriminatorByLeadingTrackPtCut)[tau] == 1)
       ) filter = true;
	else LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::filterInput:Info: Bad tau discriminator (isolation = "<<(int)(*thePFTauDiscriminatorByIsolation)[tau]<<", minTrackPt = "<<(int)(*thePFTauDiscriminatorByLeadingTrackPtCut)[tau]<<", electron veto = "<<(int)(*thePFTauDiscriminatorAgainstElectrons)[tau]<<", muon veto = "<<(int)(*thePFTauDiscriminatorAgainstMuons)[tau]<<"). Skip tauCand.";
	
	return filter;
}
reco::TrackRefVector InputTrackSelector::getPFTauDaughters(reco::PFTauRef &PFTau){
	reco::TrackRefVector trkVct;
	const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands();//cand in signal cone 
	//isolationPFChargedHadrCands stores tracks in isol/veto cone
	for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
		LogTrace("InputTrackSelector")<<"evt "<<iEvent_->id().event()<<" InputTrackSelector::getPFTauDaughters: PFTau daughter pt "<<iter->get()->pt()<<", eta "<<iter->get()->eta()<<", vtx("<<iter->get()->vx()<<","<<iter->get()->vy()<<","<<iter->get()->vz()<<")";
		if(iter->get()->trackRef().isNonnull()) trkVct.push_back( (*iter)->trackRef() );
	}
	
	return trkVct;
}


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);
