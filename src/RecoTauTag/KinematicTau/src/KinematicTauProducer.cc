#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
iConfig_(iConfig),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
usedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "usedTauCandidates" ) ),
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) )
{
	produces<int>("kinTauCreatorFlag");//0=invalid, 1=valid
	produces<InputTauCollection>("usedTauRefs");//needed to fill in unfit KinematicParticle later on
	produces<KinematicCollection>("SelectedKinematicParticles");
	produces<SelectedKinematicParticleCollection>("SelectedKinematicHiggs");//refit higgs from tau pairs
	produces<reco::PFCandidateCollection>("PFDaughters");

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
	
	std::auto_ptr<KinematicCollection> selected_ = std::auto_ptr<KinematicCollection >(new KinematicCollection);
	KinematicCollection & selected = * selected_;
	std::auto_ptr<InputTauCollection> PFTauRef_ = std::auto_ptr<InputTauCollection>(new InputTauCollection);
	InputTauCollection & PFTauRef = * PFTauRef_;
	std::auto_ptr<reco::PFCandidateCollection> PFDaughters_ = std::auto_ptr<reco::PFCandidateCollection>(new reco::PFCandidateCollection);
	reco::PFCandidateCollection & PFDaughters = * PFDaughters_;
	
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);
	
	tauM_ = new std::vector<RefCountedKinematicTree>;
	tauP_ = new std::vector<RefCountedKinematicTree>;
	
	bool filterValue;
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);
	if(primVtxs->size()<1) return false;
	const reco::Vertex primVtx = primVtxs->front();
	filterValue = select(selected, primVtx, PFTauRef, PFDaughters);

	iEvent_->put(PFTauRef_,"usedTauRefs");
	edm::OrphanHandle<reco::PFCandidateCollection> orphanPFCands = iEvent_->put(PFDaughters_,"PFDaughters");
//	correctReferences(selected, orphanPFCands);//has to be called before put(selected_,"SelectedKinematicParticles")!!!
	iEvent_->put(selected_,"SelectedKinematicParticles");
	
	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(filterValue) flag = 1;
	else flag = 0;
	iEvent_->put(flagPtr,"kinTauCreatorFlag");
	
	//combineHiggs crashes due to seg fault
	//	std::auto_ptr<SelectedKinematicParticleCollection> selectedHiggs_ = std::auto_ptr<SelectedKinematicParticleCollection >(new SelectedKinematicParticleCollection);
	//	SelectedKinematicParticleCollection & selectedHiggs = * selectedHiggs_;
	//	combineHiggs(primVtx, selectedHiggs);
	//	iEvent_->put(selectedHiggs_,"SelectedKinematicHiggs");
	
	
	delete tauM_;
	delete tauP_;
	
	return filterValue;
}
// ------------ method called once each job just before starting event loop  ------------
void KinematicTauProducer::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void KinematicTauProducer::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
	std::cout<<"=- KinematicTauProducer:: asked for >=1 kinTaus per event. efficiency = "<<ratio<<" ("<<cntFound_<<"/"<<cnt_<<")"<<std::endl;
}

bool KinematicTauProducer::select(KinematicCollection & refitParticles, const reco::Vertex & primaryVtx, InputTauCollection & PFTauRef, reco::PFCandidateCollection & PFDaughters){
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
		//		reco::PFTauRef thePFTau(usedTaus, index);//ref auf ref?
		reco::PFTauRef thePFTau = usedTaus->at(index);
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
		int fitStatus = kinTauCrtr->create(primaryVtx, input);
		if(fitStatus==1) cntValid++;
		reco::PFTau refitPFTau = kinTauCrtr->getPFTau();
		std::vector<math::XYZTLorentzVector> refitDaughters = kinTauCrtr->getRefittedChargedHadrons();
		RefCountedKinematicTree refitTree = kinTauCrtr->getKinematicTree();
		std::vector<reco::TrackRef> usedTracks = kinTauCrtr->getSelectedTracks();
//		int iterations = kinTauCrtr->iterations();
//		float csum = kinTauCrtr->csum();
	}
	if(cntValid>0) cntFound_++;//found at least one refit tau
	
//	int cntDoubleTaus = 0;
//	for (std::vector<bool>::iterator iter=ambiguity.begin(); iter!=ambiguity.end(); ++iter) {
//		if(*iter) cntDoubleTaus++;
//	}
//	if(refitParticles.size() - cntDoubleTaus >= 2){
//		cntFound_++;
//		fullyDetermined = true;
//		if(verbosity_>=2) printf("evt %d KinematicTauProducer::select: %d kinematic tau(s) reconstructed with %d ambiguities.\n", iEvent_->id().event(), refitParticles.size(), cntDoubleTaus);
//	}else printf("evt %d KinematicTauProducer::select:Warning: only %d kinematic tau(s) reconstructed with %d ambiguities. Skip Evt.\n", iEvent_->id().event(), refitParticles.size(), cntDoubleTaus);
	
	return fullyDetermined;
}

/*int ThreeProngTauCreator::addRefittedParticles(const int &ambiguityCnt, RefCountedKinematicTree tree, KinematicConstrainedVertexFitter* kcvFitter, KinematicCollection &refitParticles, reco::PFTauRef &tauRef, const reco::Vertex &primVtx){
	try{
		tree->movePointerToTheTop();
	}catch(VertexException){
		std::cout<<"KinematicTree::movePointerToTheTop; tree-> is empty! -- Event skipped."<<std::endl;
		return false;
	}
	std::vector<SelectedKinematicParticle> refitTauDecay;
	if(verbosity_>=2){
		double pt = sqrt(pow(tree->currentParticle()->currentState().kinematicParameters().momentum().x(),2)+pow(tree->currentParticle()->currentState().kinematicParameters().momentum().y(),2));
		printf("ThreeProngTauCreator::added tau: pt %f.\n", pt);
	}
	//	RefCountedKinematicParticle kinematicParticle = tree->currentParticle();
	iterations_ = kcvFitter->getIterations();
	csum_ = kcvFitter->getSumC();
	int status = 1;
	std::string name;
	const edm::ParameterSet & fitConfig = iConfig_.getParameter<edm::ParameterSet>( "fitParameter" );
	int maxiterations = fitConfig.getParameter<int>( "maxNbrOfIterations" );
	double mincsum = fitConfig.getParameter<double>( "maxDistance" );
	reco::PFCandidateRef pfref;
	
	int motherCharge = 0;
	if(tree->currentParticle()->currentState().particleCharge() < 0){
		motherCharge = -1;
		name = std::string("tauM");
	}else{
		motherCharge = 1;
		name = std::string("tauP");
	}
	refitTauDecay.push_back( SelectedKinematicParticle(tree->currentParticle(), name, iterations_, maxiterations, csum_, mincsum, pfref, ambiguityCnt, status) );
	refitTauDecay.back().setInitialTauState(tauRef, primVtx);//initial tau state consists of used primVtx (including errors) and the pfTau parameters
	
	std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
	unsigned int i = 0;
	for (std::vector<RefCountedKinematicParticle>::iterator iter=daughters.begin(); iter!=daughters.end(); ++iter, ++i) {
		if(i<=2){
			if(motherCharge == -1) name = std::string("piM")+BasicTools::intToString(i);
			else name = std::string("piP")+BasicTools::intToString(i);
		}else{
			if(motherCharge == -1) name = std::string("nuM");
			else name = std::string("nuP");
		}
		refitTauDecay.push_back( SelectedKinematicParticle(*iter, name, iterations_, maxiterations, csum_, mincsum, pfref, ambiguityCnt, status) );
	}
	
	if(refitTauDecay.size() != 5) printf("ThreeProngTauCreator::getRefittedParticles:Saved only %i refitted particles.\n", refitTauDecay.size());
	else refitParticles.push_back(refitTauDecay);
	
	return refitTauDecay.size();
	
	//	check4VctConservation(refitTau, outs);
	//	checkPrimVtxC(refitTau,*primVtx);
	//	std::cout<<"nuMass "<<refitNu->mass()<<std::endl;
}
void ThreeProngTauCreator::correctReferences(KinematicCollection & selected, edm::OrphanHandle<reco::PFCandidateCollection> & orphanPFCands){
	unsigned index = 0;
	std::vector<reco::PFCandidateRef> newPFRefs;
	for(reco::PFCandidateCollection::const_iterator iter = orphanPFCands->begin(); iter != orphanPFCands->end(); ++iter, ++index){
		reco::PFCandidateRef ref(orphanPFCands, index);
		newPFRefs.push_back(ref);
	}
	//	std::cout<<"newPFRefs size = "<<newPFRefs.size()<<std::endl;
	//	std::cout<<"KinematicCollection size = "<<selected.size()<<std::endl;
	index = 0;
	for(KinematicCollection::iterator iter = selected.begin(); iter != selected.end(); ++iter){
		if(iter->front().ambiguity() == 2) index = index-3;//if second solution the last PFRefs are used again. (2nd sol only exists if a first one exists too)
		for(SelectedKinematicParticleCollection::iterator iter2 = iter->begin(); iter2 != iter->end(); ++iter2){
			if(iter2->name().substr(0,2) == "pi"){
				if(index >= newPFRefs.size()){
					printf("evt %d ThreeProngTauCreator::correctReferences: Bad selection size! index=%d, refs=%d\n", iEvent_->id().event(), index, newPFRefs.size());
					throw 111;
				}
				//				std::cout<<"newPFRefs index = "<<index<<", pt = ";
				//				if (!newPFRefs[index].isNull()) std::cout<<newPFRefs[index]->pt()<<std::endl;
				//				else std::cout<<" NAN"<<std::endl;
				iter2->setPFCandRef(newPFRefs[index]);
				index++;
			}
		}
	}
}
 */
//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
