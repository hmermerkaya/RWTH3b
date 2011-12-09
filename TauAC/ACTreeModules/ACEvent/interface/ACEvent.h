#ifndef ACEvent_h
#define ACEvent_h
/**
 This class contains the branches of the main tree.
 
 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../../ACDataFormats/ACGeneral/interface/ACDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACEventGlobals.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACEventInfo.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFitParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACGenDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFittedDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFittedThreeProngDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACGenParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACJet.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACPFTau.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACPileupInfo.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACTrigger.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACVertex.h"

class ACEvent {
public:
    ACEvent();
    const ACEventInfo * eventInfo() const;
    /// set event info
    ACEventInfo ** linkEventInfo();
    
    const ACEventGlobals * eventGlobals() const;
    /// set event globals
    ACEventGlobals ** linkEventGlobals();
    
    const ACTrigger * trigger() const;
    /// set trigger info
    ACTrigger ** linkTrigger();
    
    const std::vector<ACVertex *>& offlinePV() const;
    /// set vertex collection for offline primay vertices
    std::vector<ACVertex *> ** linkOfflinePV();
    
    const std::vector<ACVertex *>& reducedPV() const;
    /// set vertex collection for primary vertices reduced by the kinematic tau fit
    std::vector<ACVertex *> ** linkReducedPV();
    
    const std::vector<ACGenParticle *>& generator() const;
    /// set generator collection
    std::vector<ACGenParticle *> ** linkGenerator();
    
    const std::vector<ACGenDecay *>& genTauDecays() const;
    /// set generator tau decays
    std::vector<ACGenDecay *> ** linkGenTauDecays();
    
    const std::vector<ACParticle *>& muons() const;
    /// set muon collection
    std::vector<ACParticle *> ** linkMuons();
    
    const std::vector<ACParticle *>& electrons() const;
    /// set electron collection
    std::vector<ACParticle *> ** linkElectrons();
    
    const std::vector<ACFittedThreeProngDecay *>& tauDecays() const;
    /// set tau-decay collection
    std::vector<ACFittedThreeProngDecay *> ** linkTauDecays();
    
    /// set the collection of fitted particles of the fit. will be referenced from tauDecays_. no direct get function needed.
    std::vector<ACFitParticle *> ** linkFittedThreeProngParticles();
    
    const std::vector<ACJet *>& pfJets() const;
    /// set pf jet collection
    std::vector<ACJet *> ** linkPFJets();
    
    const std::vector<ACPFTau *>& pfTaus() const;
    /// set pf tau collection
    std::vector<ACPFTau *> ** linkPFTaus();
    
protected:
    /// main branch collections
    ACEventInfo * eventInfo_;
    ACEventGlobals * eventGlobals_;
    ACTrigger * trigger_;
    std::vector<ACVertex *> * offlinePV_, * reducedPV_;
    std::vector<ACGenParticle *> * generator_;
    std::vector<ACGenDecay *> * genTauDecays_;
    std::vector<ACParticle *> * muons_, * electrons_;
    /// result of the kinematic fit
    std::vector<ACFittedThreeProngDecay *> * tauDecays_;
    /// stores the fitted particles of the fit. will be referenced from tauDecays_. no direct get function needed.
    std::vector<ACFitParticle *> * fittedThreeProngParticles_;
    std::vector<ACJet *> * pfJets_;
    std::vector<ACPFTau *> * pfTaus_;
    
    /// test the validity a branch content
    template<class T> void checkContent(const T * address, const std::string & label) const;
};


template<class T> void ACEvent::checkContent(const T * address, const std::string & label) const {
    if (address == 0) {
        printf("ACEvent::checkContent: Error! Collection %s is not linked!\n", label.c_str());
        throw 404;
    }
}

#endif
