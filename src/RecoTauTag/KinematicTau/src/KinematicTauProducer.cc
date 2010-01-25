#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
iConfig_(iConfig),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
usedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "usedTauCandidates" ) ),
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) )
{
	produces<int>("kinTauCreatorFlag");//0=invalid, 1=valid
	produces<InputTauCollection>("usedTauRefs");//needed to fill in unfit KinematicParticle later on
	produces<reco::RecoChargedCandidateCollection>("usedTauDaughters");
	produces<SelectedKinematicDecayCollection>("SelectedKinematicParticles");
	produces<SelectedKinematicParticleCollection>("SelectedKinematicHiggs");//refit higgs from tau pairs

	//SelectedKinematicParticle = 1 refitted particle
	//vector<SelectedKinematicParticle> = all particles from 1 tau decay (i.e tau->3pi+nu) starting with mother
	//vector<vector<SelectedKinematicParticle> > = several tau decays or ambiguity
}


KinematicTauProducer::~KinematicTauProducer(){
}



// ------------ method called on each new Event  ------------
bool KinematicTauProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
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
	if(primVtxs->size()<1) return false;
	const reco::Vertex primVtx = primVtxs->front();
	bool filterValue = select(selected, primVtx, PFTauRefCollection, daughterCollection);

	iEvent_->put(PFTauRefCollection_,"usedTauRefs");
	edm::OrphanHandle<reco::RecoChargedCandidateCollection> orphanCands = iEvent_->put(daughterCollection_,"usedTauDaughters");
	correctReferences(selected, orphanCands);//has to be called before put(selected_,"SelectedKinematicParticles")!!!
	iEvent_->put(selected_,"SelectedKinematicParticles");
	
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
	std::cout<<"=- KinematicTauProducer:: asked for >=1 kinTaus per event. efficiency = "<<ratio<<" ("<<cntFound_<<"/"<<cnt_<<")"<<std::endl;
}

bool KinematicTauProducer::select(SelectedKinematicDecayCollection & refitDecays, const reco::Vertex & primaryVtx, InputTauCollection & PFTauRefCollection, reco::RecoChargedCandidateCollection & daughterCollection){
	bool fullyDetermined = false;
	std::vector<bool> ambiguity;
	
	edm::Handle<InputTrackCollection> inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
	edm::Handle<InputTauCollection> usedTaus;
	iEvent_->getByLabel(usedTauCandidatesTag_, usedTaus);
	if(inputCollection->size() != usedTaus->size()){
		std::cout<<"KinematicTauProducer::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<usedTauCandidatesTag_.label()<<"("<<usedTaus->size()<<")"<<std::endl;
		return false;
	}
	unsigned int index = 0;
	TransientTrackBuilder trkBuilder = *transTrackBuilder_;
	KinematicTauCreator *kinTauCrtr = new ThreeProngTauCreator(trkBuilder, iConfig_);
	unsigned int cntValid = 0;
	for(InputTrackCollection::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks, ++index) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
		int fitStatus = kinTauCrtr->create(primaryVtx, input);
		if(fitStatus==1){
			cntValid++;
			//		reco::PFTau refitPFTau = kinTauCrtr->getPFTau();
			//		std::vector<math::XYZTLorentzVector> refitDaughters = kinTauCrtr->getRefittedChargedHadrons();
			std::vector<reco::TrackRef> usedTracks = kinTauCrtr->getSelectedTracks();
			saveSelectedTracks(usedTracks, daughterCollection);
			saveKinParticles(kinTauCrtr, refitDecays, primaryVtx);
			//save tau ref
			reco::PFTauRef tauRef = usedTaus->at(index);
			PFTauRefCollection.push_back(tauRef);
		}
	}
	if(cntValid>0){
		fullyDetermined = true;
		cntFound_++;//found at least one refit tau
	}
	
//	int cntDoubleTaus = 0;
//	for (std::vector<bool>::iterator iter=ambiguity.begin(); iter!=ambiguity.end(); ++iter) {
//		if(*iter) cntDoubleTaus++;
//	}
//	if(refitDecays.size() - cntDoubleTaus >= 2){
//		cntFound_++;
//		fullyDetermined = true;
//		if(verbosity_>=2) printf("evt %d KinematicTauProducer::select: %d kinematic tau(s) reconstructed with %d ambiguities.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
//	}else printf("evt %d KinematicTauProducer::select:Warning: only %d kinematic tau(s) reconstructed with %d ambiguities. Skip Evt.\n", iEvent_->id().event(), refitDecays.size(), cntDoubleTaus);
	
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
//	if(verbosity_>=2){
//		double pt = sqrt(pow(tree->currentParticle()->currentState().kinematicParameters().momentum().x(),2)+pow(tree->currentParticle()->currentState().kinematicParameters().momentum().y(),2));
//		printf("KinematicTauProducer::added tau: pt %f.\n", pt);
//	}
	//	RefCountedKinematicParticle kinematicParticle = tree->currentParticle();
	int iterations = kcvFitter->getNit();
	float csum = kcvFitter->getCSum();
	
	int status = 1;
	std::string name;
	int ambiguityCnt = -1;
	int maxiterations = iConfig_.getParameter<int>( "maxNbrOfIterations" );
	double mincsum = iConfig_.getParameter<double>( "maxDistance" );
	reco::RecoChargedCandidateRef emptyCandRef;//pions will be filled later on by correctReferences()
		
	int motherCharge = 0;
	if(tree->currentParticle()->currentState().particleCharge() < 0){
		motherCharge = -1;
		name = std::string("tauM");
	}else{
		motherCharge = 1;
		name = std::string("tauP");
	}
	
	SelectedKinematicParticleCollection refitTauDecay;
	refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	std::vector<reco::TrackRef> tpliter = kinTauCrtr->getSelectedTracks();
	TLorentzVector initialTauMomentum;
	for ( std::vector<reco::TrackRef>::const_iterator iter2 = tpliter.begin(); iter2 != tpliter.end(); ++iter2 ) {
		TLorentzVector p4tmp;
		p4tmp.SetVectM(TVector3((*iter2)->px(), (*iter2)->py(), (*iter2)->pz()), 0.1395702);
		initialTauMomentum += p4tmp; 
	} 
	refitTauDecay.back().setInitialTauState(initialTauMomentum, primVtx);//initial tau state consists of used primVtx (including errors) and the pfTau parameters
	
	std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
	unsigned int i = 0;
	for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter, ++i) {
		if(i<=2){
			if(motherCharge == -1) name = std::string("piM")+intToString(i);
			else name = std::string("piP")+intToString(i);
		}else{
			if(motherCharge == -1) name = std::string("nuM");
			else name = std::string("nuP");
		}
		refitTauDecay.push_back( SelectedKinematicParticle(*iter, name, iterations, maxiterations, csum, mincsum, emptyCandRef, ambiguityCnt, status) );
	}
	
	if(refitTauDecay.size() != 5) printf("KinematicTauProducer::saveSelectedTracks:Saved only %i refitted particles.\n", refitTauDecay.size());
	else refitDecays.push_back(refitTauDecay);

	return refitTauDecay.size();
	
	//	check4VctConservation(refitTau, outs);
	//	checkPrimVtxC(refitTau,*primVtx);
	//	std::cout<<"nuMass "<<refitNu->mass()<<std::endl;
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
		decay->chargedDaughters(daughters);// = decay->chargedDaughters();//needs to be copied!
//		if(decay->front().ambiguity() == 2) index = index-3;//if second solution the last PFRefs are used again. (2nd sol only exists if a first one exists too)
		for(std::vector<SelectedKinematicParticle*>::iterator particle = daughters.begin(); particle != daughters.end(); ++particle){
			if(index >= newRefs.size()){
				printf("evt %d KinematicTauProducer::correctReferences: Bad selection size! index=%d, refs=%d\n", iEvent_->id().event(), index, newRefs.size());
				throw 111;
			}
//			std::cout<<"newRefs index = "<<index<<", pt = ";
//			if (!newRefs[index].isNull()) std::cout<<newRefs[index]->pt()<<std::endl;
//			else std::cout<<" NAN"<<std::endl;
			
			if(*particle!=NULL){
//				printf("particle name: %s, pt: %f\n", (*particle)->name().c_str(), (*particle)->p4().Pt());
				(*particle)->setCandRef(newRefs[index]);
			}else std::cout<<"Reference not modified!!!(at index="<<index<<")"<<std::endl;
			index++;
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
