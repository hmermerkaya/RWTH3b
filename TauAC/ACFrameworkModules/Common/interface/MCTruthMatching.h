#ifndef MCTruthMatching_h
#define MCTruthMatching_h
/**
 Test for MCTruthMatching.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"
#include "ACDataFormats/ACGeneral/interface/ACRef.h"

class MCTruthMatching {
public:
    MCTruthMatching(edm::InputTag mapTag);
    virtual ~MCTruthMatching() {};

    /// trigger e new event. cleans up counters and references.
    void nextEvt();
    /// print a summary
    void printOutro();
    /// store the association between an offline ref and the corresponding framework ref from which it was created from.
    void logConversion(const reco::GenParticleRef& onlineRef, const ACGenParticleRef & offlineRef);
    /// find TrackRef from SelectedKinematicParticle in Association and check its generator matching.
    const ACGenParticleRef * getMatching(const edm::Event& evt, const SelectedKinematicParticle * kinparticle);
    /// retrieve a reference to the converted offline object from a framework reference
    const ACGenParticleRef * getConverted(const reco::GenParticleRef & onlineRef) const;
    /// increment counter if getMatching is not used
    void increment(bool matched);

private:
    edm::InputTag mapTag_;
    edm::Handle<reco::GenParticleMatch> map_;
    /// association between an offline ref and the corresponding framework ref from which it was created from.
    std::map<reco::GenParticleRef, ACGenParticleRef> conversionLog_;
    /// store framework reference to generator particles, which where already matched to a particle in this event
    std::vector<reco::GenParticleRef> assignedRefsPerEvt_;

    unsigned int cntMatched_, cntDuplicateMatched_;
    unsigned int passedParticles_;

    /// test whether the same generator particle was assigned more than once
    bool testDuplicate(const reco::GenParticleRef & ref);

};

#endif
