#include "../interface/InputTrackSelector.h"

InputTrackSelector::InputTrackSelector(const edm::ParameterSet& iConfig):
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauCandidates" ) )
{
	produces<int>("inputTracksFlag");//0=invalid, 1=valid
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
	iEvent_->put(flagPtr,"inputTracksFlag");
	
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
	printf("=- InputTrackSelector:: found at least 1 tau candidate per event. efficiency = %f (%i/%i)\n", ratio, cntFound_, cnt_);
}
bool InputTrackSelector::select(InputTrackCollection & selected, InputTauCollection & PFTauRef){
	bool found = false;
	
	edm::Handle<reco::PFTauCollection> inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
	for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < inputCollection->size(); iPFTau++) {
		reco::PFTauRef thePFTau(inputCollection, iPFTau);
		if(!filterInput(thePFTau)) continue;//move into external module?
		
		reco::TrackRefVector tauDaughters = getPFTauDaughters(thePFTau);
		if(tauDaughters.size()>=3){//move this number into config file?
			selected.push_back(tauDaughters);
			
			PFTauRef.push_back(thePFTau);
		}else LogTrace("KinematicTauCreator")<<"InputTrackSelector::select: only "<<tauDaughters.size()<<" tau daughter(s) found. Skip tau candidate.";
	}
	
	if(selected.size()>0){
		cntFound_++;
		found = true;
		LogTrace("KinematicTauCreator")<<"InputTrackSelector::select: "<<selected.size()<<" tau candidate(s) reconstructed.";
	}else LogTrace("KinematicTauCreator")<<"InputTrackSelector::select:Warning: only "<<selected.size()<<" tau candidate(s) reconstructed. Skip Evt.";

	return found;
}
bool InputTrackSelector::filterInput(reco::PFTauRef &tau){//use seperate filter module?
	bool filter = false;
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorByIsolation;
	iEvent_->getByLabel("fixedConeHighEffPFTauDiscriminationByTrackIsolation",thePFTauDiscriminatorByIsolation); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorByLeadingTrackPtCut; 
	iEvent_->getByLabel("fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut",thePFTauDiscriminatorByLeadingTrackPtCut); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorAgainstElectrons; 
	iEvent_->getByLabel("fixedConeHighEffPFTauDiscriminationAgainstElectron",thePFTauDiscriminatorAgainstElectrons); 
	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminatorAgainstMuons; 
	iEvent_->getByLabel("fixedConeHighEffPFTauDiscriminationAgainstMuon",thePFTauDiscriminatorAgainstMuons); 
	
//edmDumpEventContent /disk1/perchalla/data/CMSSW_3_1_2/KinTau/tau3piFromVBFH/AODSIMHLT_tau3piFromVBFH_145GeV.root | grep PFTau	
//	fixedConeHighEffPFTauDiscriminationAgainstElectron
//	fixedConeHighEffPFTauDiscriminationAgainstMuon
//	fixedConeHighEffPFTauDiscriminationByECALIsolation
//	fixedConeHighEffPFTauDiscriminationByECALIsolationUsingLeadingPion
//	fixedConeHighEffPFTauDiscriminationByIsolation
//	fixedConeHighEffPFTauDiscriminationByIsolationUsingLeadingPion
//	fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut
//	fixedConeHighEffPFTauDiscriminationByLeadingTrackFinding
//	fixedConeHighEffPFTauDiscriminationByLeadingTrackPtCut
//	fixedConeHighEffPFTauDiscriminationByTrackIsolation
//	fixedConeHighEffPFTauDiscriminationByTrackIsolationUsingLeadingPion
	
	//	edm::Handle<reco::PFTauDiscriminator> thePFTauDiscriminator!TaNC;
	//	iEvent_->getByLabel("fixedConeHighEffPFTauDiscriminationBy!TaNCfrQuarterPercent",thePFTauDiscriminator!TaNC);

//	if( (int)(*thePFTauDiscriminatorByIsolation)[tau]
//	   *(int)(*thePFTauDiscriminatorByLeadingTrackPtCut)[tau]
//	   *(int)(*thePFTauDiscriminatorAgainstElectrons)[tau]
//	   *(int)(*thePFTauDiscriminatorAgainstMuons)[tau]
//	   != 0) filter = true;
    if (
        ((*thePFTauDiscriminatorAgainstMuons)[tau] == 1 || (*thePFTauDiscriminatorAgainstElectrons)[tau] == 1)
        && ((*thePFTauDiscriminatorByIsolation)[tau] == 1 || (*thePFTauDiscriminatorByLeadingTrackPtCut)[tau] == 1)
       ) filter = true;
	else LogTrace("KinematicTauCreator")<<"InputTrackSelector::filterInput:Info: Bad tau discriminator (isolation = "<<(int)(*thePFTauDiscriminatorByIsolation)[tau]<<", minTrackPt = "<<(int)(*thePFTauDiscriminatorByLeadingTrackPtCut)[tau]<<", electron veto = "<<(int)(*thePFTauDiscriminatorAgainstElectrons)[tau]<<", muon veto = "<<(int)(*thePFTauDiscriminatorAgainstMuons)[tau]<<"). Skip tauCand.";
	
	return filter;
}
reco::TrackRefVector InputTrackSelector::getPFTauDaughters(reco::PFTauRef &PFTau){
	reco::TrackRefVector trkVct;
	const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands();//cand in signal cone 
	//isolationPFChargedHadrCands stores tracks in isol/veto cone
	for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
		LogTrace("KinematicTauCreator")<<"InputTrackSelector::getPFTauDaughters: PFTau daughter pt "<<iter->get()->pt()<<", eta "<<iter->get()->eta()<<", vtx("<<iter->get()->vx()<<","<<iter->get()->vy()<<","<<iter->get()->vz()<<")";
		if(iter->get()->trackRef().isNonnull()) trkVct.push_back( (*iter)->trackRef() );
	}
	
	return trkVct;
}


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);
