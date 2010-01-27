#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
selectedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "selectedTauCandidates" ) ),//currently not in use
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) ),
minKinTau_( iConfig.getUntrackedParameter<unsigned int>( "minKinTau", 1 ) )//filter returns true if more than minKinTau_ taus were fitted
{
	produces<int>("kinTauCreatorFlag");//0=invalid, 1=valid
	produces<InputTauCollection>("selectedTauRefs");//currently not in use
	produces<reco::RecoChargedCandidateCollection>("selectedTauDaughters");//for matching issues
	produces<SelectedKinematicDecayCollection>("SelectedKinematicDecays");
}

KinematicTauProducer::~KinematicTauProducer(){
}

// ------------ method called on each new Event  ------------
bool KinematicTauProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	bool filterValue = false;

	cnt_++;
	iEvent_ = &iEvent;
	
	std::auto_ptr<SelectedKinematicDecayCollection> selected_ = std::auto_ptr<SelectedKinematicDecayCollection >(new SelectedKinematicDecayCollection);
	SelectedKinematicDecayCollection & selected = * selected_;
	std::auto_ptr<InputTauCollection> PFTauRefCollection_ = std::auto_ptr<InputTauCollection>(new InputTauCollection);
	InputTauCollection & PFTauRefCollection = * PFTauRefCollection_;
	std::auto_ptr<reco::RecoChargedCandidateCollection> daughterCollection_ = std::auto_ptr<reco::RecoChargedCandidateCollection>(new reco::RecoChargedCandidateCollection);
	reco::RecoChargedCandidateCollection & daughterCollection = * daughterCollection_;
	
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);
	
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);

	if(primVtxs->size()>=1){
		const reco::Vertex primVtx = primVtxs->front();
		filterValue = select(selected, PFTauRefCollection, daughterCollection, primVtx);
	}
	
	iEvent_->put(PFTauRefCollection_,"selectedTauRefs");
	edm::OrphanHandle<reco::RecoChargedCandidateCollection> orphanCands = iEvent_->put(daughterCollection_,"selectedTauDaughters");
	correctReferences(selected, orphanCands);//has to be called before put(selected_,"SelectedKinematicParticles")!!!
	iEvent_->put(selected_,"SelectedKinematicDecays");
	
	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(filterValue) flag = 1;
	else flag = 0;
	iEvent_->put(flagPtr,"kinTauCreatorFlag");
	
	return filterValue;
}
void KinematicTauProducer::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}
void KinematicTauProducer::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
	std::cout<<"=- KinematicTauProducer:: asked for >= "<<minKinTau_<<" kinTaus per event. efficiency = "<<ratio<<" ("<<cntFound_<<"/"<<cnt_<<")"<<std::endl;
}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection & refitDecays, InputTauCollection & PFTauRefCollection, reco::RecoChargedCandidateCollection & daughterCollection, const reco::Vertex & primaryVtx){
	bool fullyDetermined = false;
//	std::vector<bool> ambiguity;
	
	edm::Handle<InputTrackCollection> inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
	edm::Handle<InputTauCollection> usedTaus;
	iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
	if(inputCollection->size() != usedTaus->size()){
		std::cout<<"KinematicTauProducer::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<selectedTauCandidatesTag_.label()<<"("<<usedTaus->size()<<")"<<std::endl;
		return false;
	}
	unsigned int index = 0;
	TransientTrackBuilder trkBuilder = *transTrackBuilder_;
	KinematicTauCreator *kinTauCrtr = new ThreeProngTauCreator(trkBuilder, fitParameters_);
	unsigned int cntValid = 0;
	for(InputTrackCollection::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks, ++index) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
		int fitStatus = kinTauCrtr->create(primaryVtx, input);
		if(fitStatus==1){
			cntValid++;
			
			//for PFTau issues:
			//		reco::PFTau refitPFTau = kinTauCrtr->getPFTau();
			//		std::vector<math::XYZTLorentzVector> refitDaughters = kinTauCrtr->getRefittedChargedHadrons();

			//fitting debug only:
			std::vector<reco::TrackRef> usedTracks = kinTauCrtr->getSelectedTracks();
			saveSelectedTracks(usedTracks, daughterCollection);
			saveKinParticles(kinTauCrtr, refitDecays, primaryVtx);
			//save tau ref
			reco::PFTauRef tauRef = usedTaus->at(index);
			PFTauRefCollection.push_back(tauRef);
		}
	}
	if(cntValid >= minKinTau_){
		fullyDetermined = true;
		cntFound_++;//found at least minKinTau_ refit tau
	}

	
//ambiguity issued only. still to be implemented
//	int cntDoubleTaus = 0;
//	for (std::vector<bool>::iterator iter=ambiguity.begin(); iter!=ambiguity.end(); ++iter) {
//		if(*iter) cntDoubleTaus++;
//	}
//	if(refitDecays.size() - cntDoubleTaus >= 2){
//		cntFound_++;
//		fullyDetermined = true;
//		if(verbosity_>=2) printf("evt %d KinematicTauProducer::select: %d kinematic tau(s) reconstructed with %d ambiguities.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
//	}else printf("evt %d KinematicTauProducer::select:Warning: only %d kinematic tau(s) reconstructed with %d ambiguities. Skip Evt.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
	
	delete kinTauCrtr;
	
	return fullyDetermined;
}

void KinematicTauProducer::saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection){
	//set up TrackToCandidate converter
	edm::ParameterSet pioncfg;
	pioncfg.addParameter("particleType", 211);
	converter::TrackToCandidate tk2cand(pioncfg);
	
	for (std::vector<reco::TrackRef>::const_iterator trk = usedTracks.begin(); trk != usedTracks.end(); ++trk) {
		reco::RecoChargedCandidate tmpCand;
		tk2cand.convert(*trk, tmpCand); //convert tracks to candidates
		daughterCollection.push_back(tmpCand);
	}
}
int KinematicTauProducer::saveKinParticles(KinematicTauCreator *kinTauCrtr, SelectedKinematicDecayCollection &refitDecays, const reco::Vertex &primVtx){
	RefCountedKinematicTree tree = kinTauCrtr->getKinematicTree();
	KinematicConstrainedVertexFitter *kcvFitter = kinTauCrtr->getFitter();
	try{
		tree->movePointerToTheTop();
	}catch(VertexException){
		std::cout<<"KinematicTree::movePointerToTheTop; tree is empty! -- Event skipped."<<std::endl;
		return false;
	}
	int iterations = kcvFitter->getNit();
	float csum = kcvFitter->getCSum();
	
	int status = 1;
	std::string name;
	int ambiguityCnt = -1;
	int maxiterations = fitParameters_.getParameter<int>( "maxNbrOfIterations" );
	double mincsum = fitParameters_.getParameter<double>( "maxDistance" );
	reco::RecoChargedCandidateRef emptyCandRef;//pions will be filled later on by correctReferences()
		
	if(tree->currentParticle()->currentState().particleCharge() != 0){
		name = std::string("tau");
	}else{
		printf("KinematicTauProducer::saveKinParticles: neutral tau detected. tau skipped.\n");
		return 0;
	}
	
	SelectedKinematicParticleCollection refitTauDecay;
	refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	std::vector<reco::TrackRef> selTrks = kinTauCrtr->getSelectedTracks();
	TLorentzVector initialTauMomentum;
	for(std::vector<reco::TrackRef>::const_iterator selTrk = selTrks.begin(); selTrk != selTrks.end(); ++selTrk ){
		TLorentzVector p4tmp;
		p4tmp.SetVectM(TVector3((*selTrk)->px(), (*selTrk)->py(), (*selTrk)->pz()), 0.1395702);
		initialTauMomentum += p4tmp;
	}
	refitTauDecay.back().setInitialTauState(initialTauMomentum, primVtx);//initial tau state consists of initial primVtx (including errors) and the prefit tau parameters
	
	std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
	for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter) {
		if((*iter)->currentState().particleCharge() != 0){
			name = std::string("pion");
		}else{
			name = std::string("neutrino");
		}
		refitTauDecay.push_back( SelectedKinematicParticle(*iter, name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	}
	
	if(refitTauDecay.size() != 5) printf("KinematicTauProducer::saveSelectedTracks:Saved only %i refitted particles.\n", refitTauDecay.size());
	else refitDecays.push_back(refitTauDecay);

	return refitTauDecay.size();
}
void KinematicTauProducer::correctReferences(SelectedKinematicDecayCollection & selected, edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands){
	unsigned index = 0;
	std::vector<reco::RecoChargedCandidateRef> newRefs;
//	printf("orphanCands size = %i\n", orphanCands->size());
//	printf("selected size = %i\n", selected.size());
	for(reco::RecoChargedCandidateCollection::const_iterator iter = orphanCands->begin(); iter != orphanCands->end(); ++iter, ++index){
		reco::RecoChargedCandidateRef ref(orphanCands, index);
		newRefs.push_back(ref);
	}
	index = 0;
	for(SelectedKinematicDecayCollection::iterator decay = selected.begin(); decay != selected.end(); ++decay){
		std::vector< SelectedKinematicParticle* > daughters;
		decay->chargedDaughters(daughters);
//		if(decay->front().ambiguity() == 2) index = index-3;//if second solution the last PFRefs are used again. (2nd sol only exists if a first one exists too)
		for(std::vector<SelectedKinematicParticle*>::iterator particle = daughters.begin(); particle != daughters.end(); ++particle){
			if(index >= newRefs.size()){
				printf("evt %d KinematicTauProducer::correctReferences: Bad selection size! index=%d, refs=%d\n", iEvent_->id().event(), index, newRefs.size());
				throw 111;
			}
			if(*particle!=NULL){
				(*particle)->setCandRef(newRefs[index]);
			}else std::cout<<"Reference not modified!!!(at index="<<index<<")"<<std::endl;
			index++;
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
