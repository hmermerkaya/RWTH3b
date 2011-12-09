/**
 CMSSW analyzer to store offline information.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <ACDataFormats/ACGeneral/interface/ACTools.h>

class GenSelector : public edm::EDFilter {
public:
    explicit GenSelector(const edm::ParameterSet&);
    ~GenSelector();

private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    /// decides which signal to extract from generator event.
    bool storeEvent(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
    /// searches for motherPdgID and stores its decay into tau pairs
    bool storeTauPairEvt(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef, const int & motherPdgID);
    /// searches for motherPdgID and stores its decay into a tau
    bool storeEventSingleTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
    /// store each tau followed by its stable daughters. WARNING: no mother at the beginning!
    bool storeEventAllTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
    /// searches for the decay VBF(qqH->tautau) and stores it
    bool storeEventVBFH(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef);
    /// obtain all daughters of a particle
    void findDescendents(const reco::GenParticleRef& base, std::vector<reco::GenParticleRef> & descendents, int status, int pdgId=0);
    /// store quarks (from VBF) seperately for JetAnalysis
    reco::GenParticleRefVector fillQuarkRefs(const edm::OrphanHandle<reco::GenParticleCollection> & orphanSignal);


    edm::Handle<reco::GenParticleCollection> genCandidate_;
    edm::InputTag candCollection_;
    std::string decayType_;
    int motherPdgID_;
    unsigned int cnt_, cntFound_;
    std::auto_ptr<reco::GenParticleCollection> genSignalPtr_;
    std::auto_ptr<reco::GenParticleRefVector> genSignalRefPtr_, genSignalQuarksRefPtr_;
};
