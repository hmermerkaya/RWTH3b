#include "../interface/InputTrackSelector.h"

InputTrackSelector::InputTrackSelector(const edm::ParameterSet& iConfig):
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "tauCandidates" ) ),
verbosity_( iConfig.getUntrackedParameter("verbosity", 0) )
{
	produces<int>("inputTracksFlag");//0=invalid, 1=valid
	produces<InputTrackCollection>("InputTracks");//save collection of vector<reco::CandidateRef> for each tau cand
	produces<InputTauCollection>("usedTauRefs");//needed to fill in unfit KinematicParticle later on
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

	iEvent_->put(PFTauRef_,"usedTauRefs");
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
		}else if(verbosity_>=2) printf("evt %d InputTrackSelector::select: only %d tau daughter(s) found. Skip tau candidate.\n", iEvent_->id().event(), tauDaughters.size());
	}
	
	if(selected.size()>0){
		cntFound_++;
		found = true;
		if(verbosity_>=2) printf("evt %d InputTrackSelector::select: %d tau candidate(s) reconstructed.\n", iEvent_->id().event(), selected.size());
	}else printf("evt %d InputTrackSelector::select:Warning: only %d tau candidate(s) reconstructed. Skip Evt.\n", iEvent_->id().event(), selected.size());

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

	if( (int)(*thePFTauDiscriminatorByIsolation)[tau]
	   *(int)(*thePFTauDiscriminatorByLeadingTrackPtCut)[tau]
	   *(int)(*thePFTauDiscriminatorAgainstElectrons)[tau]
	   *(int)(*thePFTauDiscriminatorAgainstMuons)[tau]
	   != 0) filter = true;
	else if(verbosity_>=2) printf("evt %d InputTrackSelector::filterInput:Info: Bad tau discriminator (isolation = %d, minTrackPt = %d, electron veto = %d, muon veto = %d). Skip tauCand.\n", iEvent_->id().event(), (int)(*thePFTauDiscriminatorByIsolation)[tau], (int)(*thePFTauDiscriminatorByLeadingTrackPtCut)[tau], (int)(*thePFTauDiscriminatorAgainstElectrons)[tau], (int)(*thePFTauDiscriminatorAgainstMuons)[tau]);
	
//	if(verbosity_>=2)
//		std::cout<<"PFDiscriminatorByIsolation value: "<<(*thePFTauDiscriminatorByIsolation)[tau] << std::endl
//		<<"PFTau Electron Discriminant value: "<<(*thePFTauDiscriminatorAgainstElectrons)[tau] << std::endl
//		<<"Pt of the PFTau "<<(*tau).pt() << std::endl;
	
	return filter;
}
reco::TrackRefVector InputTrackSelector::getPFTauDaughters(reco::PFTauRef &PFTau){
	reco::TrackRefVector trkVct;
	const reco::PFCandidateRefVector & 	cands = PFTau->signalPFChargedHadrCands();//cand in signal cone 
	//isolationPFChargedHadrCands stores tracks in isol/veto cone
	for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
		if(verbosity_>=2) std::cout<<"PFTau daughter pt "<<iter->get()->pt()<<", eta "<<iter->get()->eta()<<", vtx("<<iter->get()->vx()<<","<<iter->get()->vy()<<","<<iter->get()->vz()<<")"<<std::endl;
		if(iter->get()->trackRef().isNonnull()) trkVct.push_back( (*iter)->trackRef() );
	}
	
	return trkVct;
}


//define this as a plug-in
DEFINE_FWK_MODULE(InputTrackSelector);
