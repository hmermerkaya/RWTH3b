#include "RecoTauTag/KinematicTau/interface/KinematicTauSkim.h"

KinematicTauSkim::KinematicTauSkim(const edm::ParameterSet& iConfig):
kinTausTag_( iConfig.getParameter<edm::InputTag>( "kinematicTaus" ) ),
discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
minTau_( iConfig.getUntrackedParameter<unsigned int>( "minTau", 1 ) )//filter returns true if more than minTau_ taus are found
{
}

KinematicTauSkim::~KinematicTauSkim(){
}

bool KinematicTauSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	bool filterValue = false;
	cnt_++;

	edm::Handle<reco::PFTauCollection> tauCollection;
	iEvent.getByLabel(kinTausTag_, tauCollection);
	
	std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators;
	for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr) {
		//std::cout<< "KinematicTauSkim::filter: Discriminator name: " <<*discr<<std::endl;
		edm::Handle<reco::PFTauDiscriminator> tmpHandle;
		iEvent.getByLabel(kinTausTag_.label(), *discr, tmpHandle);
		tauDiscriminators.push_back(tmpHandle);
	}		
	
	unsigned int validTaus = 0;
	unsigned int index = 0;
	//std::cout<<"KinematicTauSkim::filter: Test "<<tauCollection->size()<<" taus"<<std::endl;
	for (reco::PFTauCollection::const_iterator tau = tauCollection->begin(); tau != tauCollection->end(); ++tau, index++) {
		bool passed = true;
		reco::PFTauRef tauRef(tauCollection, index);
		//std::cout<< "KinematicTauSkim::filter: tau at " <<index<<": (p,pt,et, vertex) ("<<tauRef->p()<<","<<tauRef->pt()<<","<<tauRef->et()<<", ("<<tauRef->vx()<<","<<tauRef->vy()<<","<<tauRef->vz()<<"))";
		for (std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator discr = tauDiscriminators.begin(); discr!=tauDiscriminators.end(); ++discr) {
			if( !((**discr)[tauRef]) ){
				//std::cout<<"\t"<<discr->provenance()->productInstanceName()<<" "<<(**discr)[tauRef]<<std::endl;
				passed = false;
				break;
			}
		}
		if(passed) validTaus++;
	}
	
	if(validTaus >= minTau_){
		cntFound_++;//found at least minTau_ refitted tau(s)
		filterValue = true;
	}
	
	return filterValue;
}

void KinematicTauSkim::beginJob(){
	cnt_ = 0;
	cntFound_ = 0;
}

void KinematicTauSkim::endJob(){
	float ratio = 0.0;
	if(cnt_!=0) ratio=(float)cntFound_/cnt_;
    edm::LogVerbatim("KinematicTauSkim")<<"--> [KinematicTauSkim] asks for >= "<<minTau_<<" tau(s) per event passing the provided discriminators. Selection efficiency: "<<cntFound_<<"/"<<cnt_<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauSkim);
