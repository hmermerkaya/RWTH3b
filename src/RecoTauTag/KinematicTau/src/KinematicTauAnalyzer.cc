/**
 Test the KinematicTau package

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"


class KinematicTauAnalyzer : public edm::EDAnalyzer {
public:
	explicit KinematicTauAnalyzer(const edm::ParameterSet&);
	~KinematicTauAnalyzer();
	
	
private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	edm::InputTag kinTausTag_;
	std::vector<std::string> discriminators_;
};

KinematicTauAnalyzer::KinematicTauAnalyzer(const edm::ParameterSet& iConfig):
kinTausTag_( iConfig.getParameter<edm::InputTag>( "kinematicTaus" ) ),
discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") )
{
}


KinematicTauAnalyzer::~KinematicTauAnalyzer(){
}

void KinematicTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	
	edm::Handle<reco::PFTauCollection> tauCollection;
	iEvent.getByLabel(kinTausTag_, tauCollection);
	std::cout<< "no. of taus = " <<tauCollection->size()<<std::endl;

	std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators;
	for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr) {
		std::cout<< "Discriminator name: " <<*discr<<std::endl;
		edm::Handle<reco::PFTauDiscriminator> tmpHandle;
		iEvent.getByLabel("KinematicTauProducer", *discr, tmpHandle);
		tauDiscriminators.push_back(tmpHandle);
	}		

	unsigned int index = 0;
	for (reco::PFTauCollection::const_iterator tau = tauCollection->begin(); tau != tauCollection->end(); ++tau, index++) {
		reco::PFTauRef tauRef(tauCollection, index);
		std::cout<< "tau at " <<index<<": (p,pt,et, vertex) ("<<tauRef->p()<<","<<tauRef->pt()<<","<<tauRef->et()<<", ("<<tauRef->vx()<<","<<tauRef->vy()<<","<<tauRef->vz()<<"))";
		for (std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator discr = tauDiscriminators.begin(); discr!=tauDiscriminators.end(); ++discr) {
			std::cout<<", "<<discr->provenance()->productInstanceName()<<" "<<(**discr)[tauRef];
		}
		std::cout<<std::endl;
	}
}

void KinematicTauAnalyzer::beginJob(){
}

void KinematicTauAnalyzer::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTauAnalyzer);
