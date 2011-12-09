/**
 CMSSW analyzer to store offline information.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */
#include <memory>
#include <TTree.h>
#include <TProcessID.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ACDataFormats/ACGeneral/interface/ACEventGlobals.h"
#include "ACDataFormats/ACGeneral/interface/ACEventInfo.h"
#include "ACDataFormats/ACGeneral/interface/ACTrigger.h"
#include "ACDataFormats/ACGeneral/interface/ACVertex.h"
#include "ACDataFormats/ACGeneral/interface/ACGenParticle.h"
#include "ACDataFormats/ACGeneral/interface/ACJet.h"
#include "ACDataFormats/ACGeneral/interface/ACPFTau.h"
#include "ACDataFormats/ACGeneral/interface/ACParticle.h"
#include "ACDataFormats/ACGeneral/interface/ACFittedThreeProngDecay.h"
#include "ACDataFormats/ACConversion/interface/ACVertexConverter.h"
#include "ACDataFormats/ACConversion/interface/ACPFTauConverter.h"
#include "ACDataFormats/ACConversion/interface/ACGenParticleConverter.h"
#include "ACDataFormats/ACConversion/interface/ACParticleConverter.h"
#include "ACDataFormats/ACConversion/interface/ACGenDecayConverter.h"
#include "ACDataFormats/ACConversion/interface/ACFittedThreeProngDecayConverter.h"
#include "ACDataFormats/ACConversion/interface/ACJetConverter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"

class FinalTreeFiller : public edm::EDAnalyzer {
public:
    explicit FinalTreeFiller(const edm::ParameterSet&);
    ~FinalTreeFiller();

private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    /// insert branches in the main tree
    void storeEvent(const edm::Event& evt);

    /// helper function to delete a pointer to a vector of pointers
    template <class T> void deleteVectorOfPointers(T * inVectorOfPointers);
    /// save way to fill a collection from a inputtag into a given handle
    template <class T> bool loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle);

    /// framework file service
    edm::Service<TFileService> fileService_;
    /// input tags
    edm::InputTag genSignalTag_, chargedTauDaughterMatchMapTag_, primVtxTag_, reducedPrimVtxTag_, pileupInfoTag_, triggerResultsTag_, muonTag_, electronTag_, kinematicTausTag_, pfMETTag_, tcMETTag_, pfJetTag_, pfTauTag_;
    std::vector<edm::InputTag> pfTauDiscriminatorTags_;
    /// event and run counters
    unsigned int evtCnt_, runCnt_;
    /// main event tree
    TTree * eventTree_;
    /// list of objects to be stored in branches
    ACEventGlobals * acEventGlobals_;
    ACEventInfo * acEventInfo_;
    ACTrigger * acTrigger_;
    std::vector<ACVertex *> * offlinePV_, * reducedPV_;
    std::vector<ACParticle *> * muons_, * electrons_;
    std::vector<ACFitParticle *> * kinematicParticles_;
    std::vector<ACFittedThreeProngDecay * > * kinematicDecays_;
    std::vector<ACGenParticle * > * genParticles_;
    std::vector<ACGenDecay * > * genTauDecays_;
    std::vector<ACJet * > * pfJets_;
    std::vector<ACPFTau * > * pfTaus_;

    HLTConfigProvider HLTCP_;
    /// true if HLT menu has changed
    bool hltChanged_;

    /// provide the association between particles from the kinematic fit and the generated particles
    MCTruthMatching * kinematicParticleMatching_;
    /// logging of the association between PFTauRefs and ACPFTauRefs
    PFTauMatching * pfTauMatching_;
};
