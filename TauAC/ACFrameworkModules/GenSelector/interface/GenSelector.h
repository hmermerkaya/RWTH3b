/**
 CMSSW analyzer searching for generator decays.

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
    /// store each tau followed by its stable daughters. if motherPdgID defined, only taus from this mother are kept. taus from different mother particles are ignored. in case of VBF Higgs, also the quarks are stored.
    bool storeEventAllTau(edm::Event& iEvent, reco::GenParticleCollection & collection, reco::GenParticleRefVector & collectionRef, int motherPdgID = 0);
    /// test whether the decay chain contains a particle with given status and pdgId
    bool hasDescendent(const reco::GenParticleRef& base, int status, int pdgId=0) const;
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
