// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      KinematicTauAdvancedProducer
// 
/**
 Entire application of KinematicTauCreator for DEBUG issues only! Use KinematicTauProducer instead.
 Part of the KinematicTau package.

 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"//own class of tauGroup to store fitted particles in event stream
#include "CommonTools/RecoAlgos/src/TrackToCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

class KinematicTauAdvancedProducer : public edm::EDFilter {
public:
    explicit KinematicTauAdvancedProducer(const edm::ParameterSet&);
    ~KinematicTauAdvancedProducer();
    
private:
    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    bool select(SelectedKinematicDecayCollection & refitDecays, reco::PFTauRefVector & PFTauRefCollection, reco::RecoChargedCandidateCollection & daughterCollection, const reco::VertexRef & primaryVtx);
    void saveSelectedTracks(const std::vector<reco::TrackRef> & usedTracks, reco::RecoChargedCandidateCollection & daughterCollection);
    int saveKinParticles(const KinematicTauCreator * kinTauCrtr, SelectedKinematicDecayCollection &refitDecays, const reco::PFTauRef & tauRef, const reco::VertexRef & primaryVtxRef);
    void correctReferences(SelectedKinematicDecayCollection & selected, const edm::OrphanHandle<reco::RecoChargedCandidateCollection> & orphanCands);
    void storePFTauDiscriminators(const reco::PFTauRef & tauRef, std::map<std::string, bool> & tauDiscriminators);
    /// store quality discriminators (depending on the decay mode)
    void setMissingQualityCriteria(SelectedKinematicDecay & decay, const KinematicTauCreator * kinTauCrtr, const reco::PFTauRef & tauRef, const reco::VertexRef & primaryVtxRef);
    
    const edm::ParameterSet fitParameters_;
    edm::Event * iEvent_;
    edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
    
    edm::InputTag primVtx_, selectedTauCandidatesTag_, inputCollectionTag_;
    std::vector<std::string> discriminators_;
    unsigned int minKinTau_, cnt_, cntFound_;
    
    static std::string intToString(int f){
        char s[32]; sprintf(s, "%d", f);
        return std::string(s);
    }    
};
