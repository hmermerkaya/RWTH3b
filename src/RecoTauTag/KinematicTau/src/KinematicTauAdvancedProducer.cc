#include "RecoTauTag/KinematicTau/interface/KinematicTauAdvancedProducer.h"

KinematicTauAdvancedProducer::KinematicTauAdvancedProducer(const edm::ParameterSet& iConfig):
fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
selectedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "selectedTauCandidates" ) ),//only used to save number of tracks in signal cone of PFTau candidate
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) ),
discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),//store the pftau discriminators
minKinTau_( iConfig.getUntrackedParameter<unsigned int>( "minKinTau", 1 ) )//filter returns true if more than minKinTau_ taus were fitted
{
	produces<int>("flag");//0=invalid, 1=valid
	produces<reco::PFTauRefVector>("selectedTauRefs");//currently not in use
	produces<reco::RecoChargedCandidateCollection>("selectedTauDaughters");//for matching issues
	produces<SelectedKinematicDecayCollection>("SelectedKinematicDecays");
}

KinematicTauAdvancedProducer::~KinematicTauAdvancedProducer(){
}

// ------------ method called on each new Event  ------------
bool KinematicTauAdvancedProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	bool filterValue = false;

	cnt_++;
	iEvent_ = &iEvent;
	
	std::auto_ptr<SelectedKinematicDecayCollection> selected_ = std::auto_ptr<SelectedKinematicDecayCollection >(new SelectedKinematicDecayCollection);
	SelectedKinematicDecayCollection & selected = * selected_;
	std::auto_ptr<reco::PFTauRefVector> PFTauRefCollection_ = std::auto_ptr<reco::PFTauRefVector>(new reco::PFTauRefVector);
	reco::PFTauRefVector & PFTauRefCollection = * PFTauRefCollection_;
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
	iEvent_->put(flagPtr,"flag");
	
	return filterValue;
}
void KinematicTauAdvancedProducer::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}
void KinematicTauAdvancedProducer::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
    edm::LogVerbatim("KinematicTauAdvancedProducer")<<"--> [KinematicTauAdvancedProducer] asks for >= "<<minKinTau_<<" kinTaus per event. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

bool KinematicTauAdvancedProducer::select(SelectedKinematicDecayCollection & refitDecays, reco::PFTauRefVector & PFTauRefCollection, reco::RecoChargedCandidateCollection & daughterCollection, const reco::Vertex & primaryVtx){
	bool fullyDetermined = false;
//	std::vector<bool> ambiguity;
	
	edm::Handle<std::vector<reco::TrackRefVector> > inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
	edm::Handle<reco::PFTauRefVector> usedTaus;
	iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
	if(inputCollection->size() != usedTaus->size()){
        edm::LogError("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<selectedTauCandidatesTag_.label()<<"("<<usedTaus->size()<<")";
		return false;
	}
	unsigned int index = 0;
	TransientTrackBuilder trkBuilder = *transTrackBuilder_;
	KinematicTauCreator *kinTauCrtr = new ThreeProngTauCreator(trkBuilder, fitParameters_);
	unsigned int cntValid = 0;
	for(std::vector<reco::TrackRefVector>::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks, ++index) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
		int fitStatus = kinTauCrtr->create(primaryVtx, input);
		if(fitStatus==1){
			cntValid++;
			
			//for PFTau issues:
			//		reco::PFTau refitPFTau = kinTauCrtr->getPFTau();
			//		std::vector<math::XYZTLorentzVector> refitDaughters = kinTauCrtr->getRefittedChargedDaughters();

			//fitting debug only:
			std::vector<reco::TrackRef> usedTracks = kinTauCrtr->getSelectedTracks();
			reco::PFTauRef tauRef = usedTaus->at(index);
			saveSelectedTracks(usedTracks, daughterCollection);
			saveKinParticles(kinTauCrtr, refitDecays, tauRef);
			//save tau ref
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
//		if(verbosity_>=2) printf("evt %d KinematicTauAdvancedProducer::select: %d kinematic tau(s) reconstructed with %d ambiguities.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
//	}else printf("evt %d KinematicTauAdvancedProducer::select:Warning: only %d kinematic tau(s) reconstructed with %d ambiguities. Skip Evt.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
	
	delete kinTauCrtr;
	
	return fullyDetermined;
}

void KinematicTauAdvancedProducer::saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection){
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
int KinematicTauAdvancedProducer::saveKinParticles(KinematicTauCreator *kinTauCrtr, SelectedKinematicDecayCollection &refitDecays, const reco::PFTauRef & tauRef){
	RefCountedKinematicTree tree = kinTauCrtr->getKinematicTree();
	KinematicConstrainedVertexFitter *kcvFitter = kinTauCrtr->getFitter();
	try{
		tree->movePointerToTheTop();
	}catch(VertexException){
		edm::LogWarning("KinematicTauAdvancedProducer")<<"KinematicTree::movePointerToTheTop; tree is empty! -- Event skipped.";
		return false;
	}
	int iterations = kcvFitter->getNit();
	float csum = kcvFitter->getCSum();
	
	int status = 1;
	std::string name;
	int ambiguityCnt = -1;
	int maxiterations = fitParameters_.getParameter<int>( "maxNbrOfIterations" );
	double mincsum = fitParameters_.getParameter<double>( "maxDelta" );
	reco::RecoChargedCandidateRef emptyCandRef;//pions will be filled later on by correctReferences()
		
	if(tree->currentParticle()->currentState().particleCharge() != 0){
		name = std::string("tau");
	}else{
		LogTrace("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::saveKinParticles: neutral tau detected. tau skipped.";
		return 0;
	}
	
	SelectedKinematicParticleCollection refitTauDecay;
	refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	refitTauDecay.back().setInitialState(TLorentzVector(tauRef->px(), tauRef->py(), tauRef->pz(), tauRef->energy()), kinTauCrtr->getModifiedPrimaryVertex());//initial tau state consists of rotated primVtx (including initial errors) and the pftau parameters
	
	std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
	for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter) {
		if((*iter)->currentState().particleCharge() != 0){
			name = std::string("pion");
		}else{
			name = std::string("neutrino");
		}
		refitTauDecay.push_back( SelectedKinematicParticle(*iter, name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	}
	
	std::map<std::string, bool> tauDiscriminators;
	storePFTauDiscriminators(tauRef, tauDiscriminators);
	
	if(refitTauDecay.size() != 5) LogTrace("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::saveSelectedTracks:Saved only "<<refitTauDecay.size()<<" refitted particles.";
	else{
		refitDecays.push_back( SelectedKinematicDecay(refitTauDecay, tauRef->signalPFChargedHadrCands().size(), tauRef->signalPFNeutrHadrCands().size(), tauDiscriminators) );
		refitDecays.back().addPFTauRef(tauRef.index());
	}
	
	return refitTauDecay.size();
}
void KinematicTauAdvancedProducer::correctReferences(SelectedKinematicDecayCollection & selected, edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands){
	unsigned index = 0;
	std::vector<reco::RecoChargedCandidateRef> newRefs;
	for(reco::RecoChargedCandidateCollection::const_iterator iter = orphanCands->begin(); iter != orphanCands->end(); ++iter, ++index){
		reco::RecoChargedCandidateRef ref(orphanCands, index);
		newRefs.push_back(ref);
	}
	index = 0;
	for(SelectedKinematicDecayCollection::iterator decay = selected.begin(); decay != selected.end(); ++decay){
		std::vector< SelectedKinematicParticle* > daughters;
		decay->modifiableChargedDaughters(daughters);
//		if(decay->front().ambiguity() == 2) index = index-3;//if second solution the last PFRefs are used again. (2nd sol only exists if a first one exists too)
		for(std::vector<SelectedKinematicParticle*>::iterator particle = daughters.begin(); particle != daughters.end(); ++particle){
			if(index >= newRefs.size()){
                edm::LogError("KinematicTauAdvancedProducer")<<"evt "<<iEvent_->id().event()<<" KinematicTauAdvancedProducer::correctReferences: Bad selection size! index="<<index<<", refs="<<newRefs.size();
				throw 111;
			}
			if(*particle!=NULL){
				(*particle)->setCandRef(newRefs[index]);
			}else edm::LogError("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::correctReferences: Reference not modified!!!(at index="<<index<<")";
			index++;
		}
	}
}

void KinematicTauAdvancedProducer::storePFTauDiscriminators(const reco::PFTauRef & tauRef, std::map<std::string, bool> & tauDiscriminators){
	//store pftau discriminators for each SelectedKinematicDecay
	if(tauDiscriminators.size()!=0){
		edm::LogWarning("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::storePFTauDiscriminators:Warning! Please provide a clean map and not of size "<<tauDiscriminators.size();
		tauDiscriminators.clear();
	}
	
	for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr) {
		edm::Handle<reco::PFTauDiscriminator> tmpHandle;
		iEvent_->getByLabel(*discr, tmpHandle);
		
		if(tauDiscriminators.find(*discr)!=tauDiscriminators.end()){
			edm::LogWarning("KinematicTauAdvancedProducer")<<"KinematicTauAdvancedProducer::storePFTauDiscriminators:Warning! Duplicate found in list of discriminators ("<<discr->c_str()<<"). Please correct this in your cfi.";
			continue;
		}
		tauDiscriminators.insert( std::make_pair( *discr, (*tmpHandle)[tauRef] ) );
	}
}
//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauAdvancedProducer);
