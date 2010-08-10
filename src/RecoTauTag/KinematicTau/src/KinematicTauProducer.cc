#include "RecoTauTag/KinematicTau/interface/KinematicTauProducer.h"

KinematicTauProducer::KinematicTauProducer(const edm::ParameterSet& iConfig):
fitParameters_( iConfig.getParameter<edm::ParameterSet>( "fitParameters" ) ),
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
selectedTauCandidatesTag_( iConfig.getParameter<edm::InputTag>( "selectedTauCandidates" ) ),//only used to save number of tracks in signal cone of PFTau candidate
inputCollectionTag_( iConfig.getParameter<edm::InputTag>( "inputTracks" ) )
{
	produces<reco::PFTauCollection>();
	produces<reco::PFTauDiscriminator>("PFRecoTauDiscriminationByKinematicFit");
	produces<reco::PFTauDiscriminator>("PFRecoTauDiscriminationByKinematicFitQuality");
}

KinematicTauProducer::~KinematicTauProducer(){
}

// ------------ method called on each new Event  ------------
bool KinematicTauProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	bool filterValue = false;

	std::auto_ptr<reco::PFTauCollection> selected_ = std::auto_ptr<reco::PFTauCollection >(new reco::PFTauCollection);
	reco::PFTauCollection & selected = * selected_;
	
	iEvent_ = &iEvent;
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);

	//make copy of initial tau collection with unmodified 4vects
	edm::Handle<InputTauCollection> usedTaus;
	iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
	selected.insert(selected.begin(), usedTaus->product()->begin(), usedTaus->product()->end());

	std::map<int, std::vector<bool> > discrimValues;
	if(primVtxs->size()>=1){
		iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);
		const reco::Vertex primVtx = primVtxs->front();

		filterValue = select(selected, discrimValues, primVtx);
	}
	
	edm::OrphanHandle<reco::PFTauCollection> orphanTaus = iEvent_->put(selected_);
	discriminate(orphanTaus, discrimValues);
	
	return filterValue;
}

void KinematicTauProducer::beginJob(){
}
void KinematicTauProducer::endJob(){
}

bool KinematicTauProducer::select(reco::PFTauCollection & selected, std::map<int, std::vector<bool> > & discrimValues, const reco::Vertex & primaryVtx){
	bool success = false;
	discrimValues.clear();
	
	edm::Handle<InputTrackCollection> inputCollection;
	iEvent_->getByLabel(inputCollectionTag_, inputCollection);
	edm::Handle<InputTauCollection> usedTaus;
	iEvent_->getByLabel(selectedTauCandidatesTag_, usedTaus);
	if(inputCollection->size() != usedTaus->size()){
		edm::LogError("KinematicTauProducer")<<"KinematicTauProducer::select: Bad input collections. Size mismatch between "<<inputCollectionTag_.label()<<"("<<inputCollection->size()<<") and "<<selectedTauCandidatesTag_.label()<<"("<<usedTaus->size()<<")";
		return false;
	}
	unsigned int index = 0;
	TransientTrackBuilder trkBuilder = *transTrackBuilder_;
	KinematicTauCreator *kinTauCrtr = new ThreeProngTauCreator(trkBuilder, fitParameters_);
	for(InputTrackCollection::const_iterator tracks = inputCollection->begin(); tracks != inputCollection->end(); ++tracks, ++index) {
		std::vector<reco::TrackRef> input;
		for(reco::TrackRefVector::iterator trk = tracks->begin(); trk!=tracks->end(); ++trk) input.push_back(*trk);
		int fitStatus = kinTauCrtr->create(primaryVtx, input);

		//save std::map<int, std::vector<bool> >, stores association between position in tauCollection and one bool for every discriminator
		reco::PFTauRef tauRef = usedTaus->at(index);
		discrimValues.insert(std::make_pair(tauRef.index(), std::vector<bool>()));
		discrimValues.find(tauRef.index())->second.push_back(fitStatus);
		discrimValues.find(tauRef.index())->second.push_back(dicriminatorByKinematicFitQuality(kinTauCrtr, fitStatus, tauRef, primaryVtx));
		
		if(fitStatus==1){
			success = true;
			//modify tau in selected list
			reco::PFTau refitPFTau = kinTauCrtr->getPFTau(); //this is only the visible part of the tau momentum!
			selected.at(tauRef.index()).setP4(refitPFTau.p4());
			selected.at(tauRef.index()).setVertex(refitPFTau.vertex()); //this is the rotated primary vertex
		}
	}
	
	delete kinTauCrtr;
	
	return success; //at least one tau was fitted
}
bool KinematicTauProducer::dicriminatorByKinematicFitQuality(const KinematicTauCreator *kinTauCrtr, const int & fitStatus, const reco::PFTauRef & tauRef, const reco::Vertex & primaryVtx){
	//combine a discriminator of important quality cuts
	
	//test if fit could create the final decay tree
	if(!fitStatus) return false;
	reco::PFTau refitPFTau = kinTauCrtr->getPFTau();//only the visible part!
	std::vector<math::XYZTLorentzVector> chargedDaughters = kinTauCrtr->getRefittedChargedDaughters();
	std::vector<math::XYZTLorentzVector> neutralDaughters = kinTauCrtr->getRefittedNeutralDaughters();
    
	//vertex separation between the modified primary vertex and the secondary vertex obtained by the fit
	reco::Vertex modifiedPV = kinTauCrtr->getModifiedPrimaryVertex();
    VertexState secVtx(kinTauCrtr->getKinematicTree()->currentDecayVertex()->position(), kinTauCrtr->getKinematicTree()->currentDecayVertex()->error());
    VertexDistance3D vtxdist;
    if ( vtxdist.distance(modifiedPV, secVtx).value() < 0.1 ) return false;
	
	//vertex significance between modified and initial primary vertex
	if ( vtxdist.distance(modifiedPV, primaryVtx).significance() > 3.0 ) return false;
	
	//WARNING!!!
	//from now one we assume a tau decay into three pions and neutrino
	//other channels need their one discriminators
	//!!!
	if(chargedDaughters.size()!=3 || neutralDaughters.size()!=1){
		edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::dicriminatorByKinematicFitQuality:WARNING!!! KinematicTauProducer assumes a tau decay into three pions and neutrino but recieved "<<chargedDaughters.size()<<" charged and "<<neutralDaughters.size()<<" neutral daughters!";
		return false;
	}
	
	//tracks in signal cone of initial pftau candidate
    if ( tauRef->signalPFChargedHadrCands().size() > 5 ) return false;

    //dR between fitted a1 and neutrino
	math::XYZTLorentzVector sum;//assume a1=sum(3pi)
	for(std::vector<math::XYZTLorentzVector>::const_iterator daughter = chargedDaughters.begin(); daughter!=chargedDaughters.end(); ++daughter){
		sum += *daughter;
	}
	TLorentzVector a1;
	a1.SetXYZM(sum.px(), sum.py(), sum.pz(), sum.M());
	TLorentzVector nu;
	nu.SetXYZM(neutralDaughters.front().px(), neutralDaughters.front().py(), neutralDaughters.front().pz(), neutralDaughters.front().M());
	if( a1.DeltaR(nu) > .09)  return false;
    
	//maximal allowed value of GJ angle
	if( thetaGJMax(a1.M(), a1.P()) > .02) return false;

	//a1 mass
	if( a1.M() < 0.6 || a1.M() > 1.8 ) return false;
	
	return true;
}
void KinematicTauProducer::discriminate(const edm::OrphanHandle<reco::PFTauCollection> & collection, const std::map<int, std::vector<bool> > & discrimValues){
	std::auto_ptr<reco::PFTauDiscriminator> discrKinFit = std::auto_ptr<reco::PFTauDiscriminator>(new reco::PFTauDiscriminator(reco::PFTauRefProd(collection)));
	std::auto_ptr<reco::PFTauDiscriminator> discrKinFitQual = std::auto_ptr<reco::PFTauDiscriminator>(new reco::PFTauDiscriminator(reco::PFTauRefProd(collection)));
		
	for(std::map<int, std::vector<bool> >::const_iterator iter = discrimValues.begin(); iter != discrimValues.end(); ++iter){
		reco::PFTauRef tauRef(collection, iter->first);
		if(iter->second.size() < 2){
			edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::discriminate: ERROR! Bad discriminator size. This tau will be skipped.";
			continue;
		}
		discrKinFit->setValue(iter->first, iter->second.at(0));
		discrKinFitQual->setValue(iter->first, iter->second.at(1));
	}
	
	iEvent_->put(discrKinFit, "PFRecoTauDiscriminationByKinematicFit");
	iEvent_->put(discrKinFitQual, "PFRecoTauDiscriminationByKinematicFitQuality");
	
}

double KinematicTauProducer::thetaGJMax(double ma1, double pa1, double Mtau){
	double argument = (-pow(ma1,2.) + pow(Mtau,2.))/(2.*Mtau*pa1);// can be negative
	//catch nan
	if(fabs(argument) >  1.0) edm::LogWarning("KinematicTauProducer")<<"KinematicTauProducer::thetaGJMax: Warning! arcsin("<<argument<<") = "<<asin(argument)<<". (pa1 "<<pa1<<", ma1 "<<ma1<<")";
	if(argument >  1.0) argument =  1.0;
	if(argument < -1.0) argument = -1.0;
	
	return asin(argument);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauProducer);
