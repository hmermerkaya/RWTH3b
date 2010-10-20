#include "../interface/GenSelectorGF.h"

GenSelectorGF::GenSelectorGF(const edm::ParameterSet& iConfig):
candCollection_( iConfig.getParameter<edm::InputTag>( "candCollection" ) ),
decay_(iConfig.getUntrackedParameter("decayType",std::string("Z3pr")))
{
	produces<int>("flag");//0=invalid, 1=valid
	
	produces<reco::GenParticleCollection>("genSignalDecay");//fixed order and size
	produces<reco::GenParticleRefVector>("genSignalDecayRef");//store refs to input list to access signal particles in matchMap
}


GenSelectorGF::~GenSelectorGF(){
}



// ------------ method called on each new Event  ------------
bool GenSelectorGF::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	cnt_++;
	
	iEvent.getByLabel( candCollection_, genCandidate_ );
	std::auto_ptr<reco::GenParticleCollection> genSignalPtr = std::auto_ptr<reco::GenParticleCollection>(new reco::GenParticleCollection);
	std::auto_ptr<reco::GenParticleRefVector> genSignalRefPtr = std::auto_ptr<reco::GenParticleRefVector >(new reco::GenParticleRefVector);
	
	reco::GenParticleCollection & genSignal = * genSignalPtr;
	reco::GenParticleRefVector & genSignalRef = * genSignalRefPtr;
	
	bool filterValue = checkGenEvt(iEvent, genSignal, genSignalRef);
	
	iEvent.put(genSignalPtr,"genSignalDecay");
	iEvent.put(genSignalRefPtr,"genSignalDecayRef");
	
	
	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(filterValue) flag = 1;
	else flag = 0;
	iEvent.put(flagPtr,"flag");
	
	return filterValue;
}

// ------------ method called once each job just before starting event loop  ------------
void GenSelectorGF::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void GenSelectorGF::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
	if(decay_!="unknown") printf("--> [GenSelectorGF] found at least two opp. charged taus. Efficiency: %d/%d = %.2f%%\n", cntFound_, cnt_, ratio*100.0);
	else printf("--> [GenSelectorGF] found at least one unknown tau decay. Efficiency: %d/%d = %.2f%%\n", cntFound_, cnt_, ratio*100.0);
}

bool GenSelectorGF::checkGenEvt(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef){
	if (decay_=="Z3pr") return checkGenEvtZ3pr(iEvent, collection, collectionRef);
	
	std::cout<<"This decay is not implemented yet."<<std::endl;
	return false;
}
bool GenSelectorGF::checkGenEvtZ3pr(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef){
	bool gen = false;
	std::vector<const reco::GenParticle* > mother;
	std::vector<const reco::GenParticle* > tauP;//charge+
	std::vector<const reco::GenParticle* > tauM;//charge-
	std::vector<const reco::GenParticle* > piP;
	std::vector<const reco::GenParticle* > piM;
	std::vector<const reco::GenParticle* > nuP;
	std::vector<const reco::GenParticle* > nuM;
	
	std::vector<reco::GenParticleRef> motherRef;
	std::vector<reco::GenParticleRef> tauPRef;//charge+
	std::vector<reco::GenParticleRef> tauMRef;//charge-
	std::vector<reco::GenParticleRef> piPRef;
	std::vector<reco::GenParticleRef> piMRef;
	std::vector<reco::GenParticleRef> nuPRef;
	std::vector<reco::GenParticleRef> nuMRef;
	
	unsigned index = 0;
	for(reco::GenParticleCollection::const_iterator candIter = genCandidate_->begin(); candIter != genCandidate_->end(); ++candIter, ++index){
		if ( abs(candIter->pdgId())==23 && candIter->numberOfDaughters()==3 ) {//Z->tautau(Z)
			mother.push_back( &*candIter );
			reco::GenParticleRef ref( genCandidate_, index );
			motherRef.push_back(ref);
		}
		if ( abs(candIter->pdgId()) == 15 && candIter->numberOfDaughters()==4) {
			const reco::Candidate *genMother = (*candIter).mother();
			if(!genMother) continue;
			genMother = ((*candIter).mother())->mother();
			if(!genMother && mother.size()>0 ) continue;
			if(genMother == mother.at(0)){
				int cntDaughterChecked = 0;
				for(unsigned int k = 0; k!=candIter->numberOfDaughters(); k++){
					if(abs(candIter->daughter(k)->pdgId())!=16 && abs(candIter->daughter(k)->pdgId())!=211) break;
					cntDaughterChecked++;
				}				
				if(cntDaughterChecked!=4) continue;
				if (candIter->pdgId()==15){
					tauM.push_back( &*candIter);
					reco::GenParticleRef ref( genCandidate_, index );
					tauMRef.push_back(ref);
				}else{
					tauP.push_back( &*candIter);
					reco::GenParticleRef ref( genCandidate_, index );
					tauPRef.push_back(ref);
				}
			}else std::cout<<"GenSelectorGF:: found other tau decaying into 3pi+nu."<<std::endl;
		}
		//		std::cout<<"taus: "<<tauM.size()+tauP.size()<<std::endl;
		if(tauM.size()==0 && tauP.size()==0) continue;
		if ( abs(candIter->pdgId()) == 16){
			const reco::Candidate *genMother = (*candIter).mother();
			if(!genMother) continue;
			if(candIter->mother()->pdgId() == 15) if(candIter->mother() == tauM.at(0)){
				nuM.push_back( &*candIter );
				reco::GenParticleRef ref( genCandidate_, index );
				nuMRef.push_back(ref);
			}
			if(candIter->mother()->pdgId() ==-15) if(candIter->mother() == tauP.at(0)){
				nuP.push_back( &*candIter );
				reco::GenParticleRef ref( genCandidate_, index );
				nuPRef.push_back(ref);
			}
		}
		if ( abs(candIter->pdgId()) == 211){
			const reco::Candidate *genMother = (*candIter).mother();
			if(!genMother) continue;
			if(candIter->mother()->pdgId() == 15) if(candIter->mother() == tauM.at(0)){
				piM.push_back( &*candIter );
				reco::GenParticleRef ref( genCandidate_, index );
				piMRef.push_back(ref);
			}
			if(candIter->mother()->pdgId() ==-15) if(candIter->mother() == tauP.at(0)){
				piP.push_back( &*candIter );
				reco::GenParticleRef ref( genCandidate_, index );
				piPRef.push_back(ref);
			}
		}
	}
	if(mother.size()!=1){
		if(mother.size()<1) std::cout<<"No Z detected.";
		if(mother.size()>1) std::cout<<"More than one Z detected.";
		std::cout<<" Skip event."<<std::endl;
		return false;
	}
	if (tauP.size()==1 && tauM.size()==1 && nuP.size()==1 && nuM.size()==1 && piP.size()==3 && piM.size()==3){
		//		std::cout<<"H->tau+tau->2x3prong found";
		cntFound_++;
		gen = true;
		
		//fixed order and size!
		collection.push_back(*(mother.at(0)));
		collection.push_back(*(tauM.at(0)));
		for(unsigned int i=0; i!=piM.size(); i++) collection.push_back(*(piM.at(i)));
		collection.push_back(*(nuM.at(0)));
		collection.push_back(*(tauP.at(0)));
		for(unsigned int i=0; i!=piP.size(); i++) collection.push_back(*(piP.at(i)));
		collection.push_back(*(nuP.at(0)));
		
		
		collectionRef.push_back(motherRef.at(0));
		collectionRef.push_back(tauMRef.at(0));
		for(unsigned int i=0; i!=piMRef.size(); i++) collectionRef.push_back(piMRef.at(i));
		collectionRef.push_back(nuMRef.at(0));
		collectionRef.push_back(tauPRef.at(0));
		for(unsigned int i=0; i!=piPRef.size(); i++) collectionRef.push_back(piPRef.at(i));
		collectionRef.push_back(nuPRef.at(0));
		//		printf("nuM = (%8.6f, %8.6f, %8.6f)\n", nuMRef.at(0)->vx(), nuMRef.at(0)->vy(), nuMRef.at(0)->vz());
		//		printf("nuP = (%8.6f, %8.6f, %8.6f)\n", nuPRef.at(0)->vx(), nuPRef.at(0)->vy(), nuPRef.at(0)->vz());
	}else{
		std::cout<<"evt"<<iEvent.id().event()<<" GenSelectorGF::checkGenEvtZ3pr: Z->tau+tau->2x3prong not found\t";
		//		std::cout<<"#of tau "<<tauCand->size()<<"#of nu "<<nuCand->size()<<"#of pi "<<piCand->size()<<std::endl;
		//	std::cout<<"evtType"<<std::endl;
		std::cout<<"sizes: nus "<<nuM.size()<<nuP.size()<<", pis "<<piM.size()<<piP.size()<<std::endl;
		std::cout<<"gen-";
		gen = false;
	}
	
	return gen;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSelectorGF);
