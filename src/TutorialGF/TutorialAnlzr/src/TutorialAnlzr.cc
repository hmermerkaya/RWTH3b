#include "../interface/TutorialAnlzr.h"

TutorialAnlzr::TutorialAnlzr(const edm::ParameterSet& iConfig):
generatorTag_( iConfig.getParameter<edm::InputTag>( "generator" ) )
{
}

TutorialAnlzr::~TutorialAnlzr(){
}

void TutorialAnlzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	edm::Handle<reco::GenParticleCollection> generatorCollection;
	iEvent.getByLabel(generatorTag_, generatorCollection);
	std::cout<< "no. of selected particles from generator = " <<generatorCollection->size()<<std::endl;

	unsigned int index = 0;
	for (reco::GenParticleCollection::const_iterator particle = generatorCollection->begin(); particle != generatorCollection->end(); ++particle, index++) {
		std::cout<< "particle at " <<index<<": (p,pt,et, vertex) ("<<particle->p()<<","<<particle->pt()<<","<<particle->et()<<", ("<<particle->vx()<<","<<particle->vy()<<","<<particle->vz()<<"))"<<std::endl;
	}
}

void TutorialAnlzr::beginJob(){
}

void TutorialAnlzr::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(TutorialAnlzr);
