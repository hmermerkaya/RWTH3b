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
#include "../../../ACDataFormats/ACGeneral/interface/ACEventWeight.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFitParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACGenDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFittedDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACFittedThreeProngDecay.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACGenParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACJet.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACPFTau.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACMuon.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACParticle.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACPileupInfo.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACTrack.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACTrigger.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACTriggerObject.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACBeamSpot.h"
#include "../../../ACDataFormats/ACGeneral/interface/ACVertex.h"

class ACEvent {
public:
    ACEvent();
    const ACEventInfo * eventInfo() const;
    /// set event info
    ACEventInfo ** linkEventInfo();

    const ACEventWeight * eventWeight() const;
    /// set event weight
    ACEventWeight ** linkEventWeight();
    
    const ACEventGlobals * eventGlobals() const;
    /// set event globals
    ACEventGlobals ** linkEventGlobals();
    
    const ACTrigger * trigger() const;
    /// set trigger info
    ACTrigger ** linkTrigger();

    const std::vector<ACTriggerObject *>& triggerObjects() const;
    /// set trigger objects
    std::vector<ACTriggerObject *> ** linkTriggerObjects();

    const ACBeamSpot * beamSpot() const;
    /// set vertex collection for offline primay vertices
    ACBeamSpot ** linkBeamSpot();

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
    
    const std::vector<ACTrack *>& tracks() const;
    /// set track collection
    std::vector<ACTrack *> ** linkTracks();

    const std::vector<ACMuon *>& muons() const;
    /// set muon collection
    std::vector<ACMuon *> ** linkMuons();
    
    const std::vector<ACParticle *>& electrons() const;
    /// set electron collection
    std::vector<ACParticle *> ** linkElectrons();
    
    const std::vector<ACFittedThreeProngDecay *>& tauDecays() const;
    /// set tau-decay collection
    std::vector<ACFittedThreeProngDecay *> ** linkTauDecays();
    
    /// set the collection of fitted particles of the fit. they will be referenced from tauDecays_. no direct get function needed.
    std::vector<ACFitParticle *> ** linkFittedThreeProngParticles();
    
    const std::vector<ACJet *>& pfJets() const;
    /// set pf jet collection
    std::vector<ACJet *> ** linkPFJets();
    
    const std::vector<ACPFTau *>& pfTaus() const;
    /// set pf tau collection
    std::vector<ACPFTau *> ** linkPFTaus();

    const std::vector<ACPileupInfo *>& pileup() const;
    /// set pileup information
    std::vector<ACPileupInfo *> ** linkPileup();

protected:
    /// branch content: basic event information
    ACEventInfo * eventInfo_;
    /// branch content: event weight
    ACEventWeight * eventWeight_;
    /// branch content: global event variables like MET
    ACEventGlobals * eventGlobals_;
    /// branch content: HLT trigger menu
    ACTrigger * trigger_;
    /// branch content: HLT trigger menu
    std::vector<ACTriggerObject *> * triggerObjects_;
    /// branch content: beam spot
    ACBeamSpot * beamSpot_;
    /// branch content: offline primary vertices
    std::vector<ACVertex *> * offlinePV_;
    /// branch content: primary vertices after removal of tracks assigned to tau daughters
    std::vector<ACVertex *> * reducedPV_;
    /// branch content: track collection
    std::vector<ACTrack *> * tracks_;
    /// branch content: muon collection
    std::vector<ACMuon *> * muons_;
    /// branch content: electron collection (gsf)
    std::vector<ACParticle *> * electrons_;
    /// branch content: particle-flow jets
    std::vector<ACJet *> * pfJets_;
    /// branch content: particle-flow taus
    std::vector<ACPFTau *> * pfTaus_;
    /// branch content: fitted decay as result of the kinematic fit. reference to assigned particles available.
    std::vector<ACFittedThreeProngDecay *> * tauDecays_;
    /// branch content: fitted particles from the kinematic tau fit. they will be referenced from tauDecays_. no direct get function needed.
    std::vector<ACFitParticle *> * fittedThreeProngParticles_;
    /// branch content: selection of generator particles
    std::vector<ACGenParticle *> * generator_;
    /// branch content: generator tau decays. reference to assigned generator particles available. 
    std::vector<ACGenDecay *> * genTauDecays_;
    /// branch content: pileup information for each bunch crossing
    std::vector<ACPileupInfo *> * pileup_;
    
    /// tests the validity a branch content
    template<class T> void checkContent(const T * address, const std::string & label) const;
};


template<class T> void ACEvent::checkContent(const T * address, const std::string & label) const {
    if (address == 0) {
        printf("ACEvent::checkContent: Error! Collection %s is not linked!\n", label.c_str());
        throw 404;
    }
}

#endif
