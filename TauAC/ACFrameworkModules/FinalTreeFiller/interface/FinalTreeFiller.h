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
#include "ACDataFormats/ACGeneral/interface/ACEventWeight.h"
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
#include "ACDataFormats/ACConversion/interface/ACPileupConverter.h"
#include "ACDataFormats/ACConversion/interface/ACTrackConverter.h"

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

#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

class FinalTreeFiller : public edm::EDAnalyzer {
public:
    explicit FinalTreeFiller(const edm::ParameterSet&);
    ~FinalTreeFiller();

private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    /// insert branches in the main tree
    void storeEvent(const edm::Event& evt);

    /// helper function to delete pointers contained in a vector pointer
    template <class T> void deleteVectorOfPointers(T * inVectorOfPointers);
    /// save way to fill a collection from a inputtag into a given handle
    template <class T> bool loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle, bool silent = false);

    /// framework file service
    edm::Service<TFileService> fileService_;
    /// input tags
    edm::InputTag genSignalTag_, genSignalRefTag_, chargedTauDaughterMatchMapTag_, primVtxTag_, reducedPrimVtxTag_, pileupInfoTag_, triggerResultsTag_, trackTag_, muonTag_, electronTag_, kinematicTausTag_, pfMETTag_, tcMETTag_, pfJetTag_, pfTauTag_;
    std::vector<edm::InputTag> pfTauDiscriminatorTags_;
    std::vector<std::string> flags_;
    std::string pileUpDistributionFileMC_, pileUpDistributionHistMC_, pileUpDistributionFileData_, pileUpDistributionHistData_;

    /// event and run counters
    unsigned int evtCnt_, runCnt_;
    /// main event tree
    TTree * eventTree_;
    
    /// branch content: basic event information
    ACEventInfo * eventInfo_;
    /// branch content: event-weight information
    ACEventWeight * eventWeight_;
    /// branch content: global event variables like MET
    ACEventGlobals * eventGlobals_;
    /// branch content: HLT trigger menu
    ACTrigger * trigger_;
    /// branch content: offline primary vertices
    std::vector<ACVertex *> * offlinePV_;
    /// branch content: primary vertices after removal of tracks assigned to tau daughters
    std::vector<ACVertex *> * reducedPV_;
    /// branch content: (filtered) tracks collection
    std::vector<ACTrack *> * tracks_;
    /// branch content: muon collection
    std::vector<ACParticle *> * muons_;
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

    HLTConfigProvider HLTCP_;
    /// true if HLT menu has changed
    bool hltChanged_;

    /// provide the association between particles from the kinematic fit and the generated particles
    MCTruthMatching * kinematicParticleMatching_;
    /// logging of the association between CMSSW PFTauRefs and ACPFTauRefs
    PFTauConversionLog * conversionLogPFTau_;
    /// logging of the association between CMSSW JetRefs and ACJetRefs
    PFJetConversionLog * conversionLogPFJet_;
    
    edm::Lumi3DReWeighting * lumiWeights_; 
};

template <class T> void FinalTreeFiller::deleteVectorOfPointers(T * inVectorOfPointers) {
    if (!inVectorOfPointers) {
        return;
    }
    typename T::iterator i;
    for (i = inVectorOfPointers->begin(); i < inVectorOfPointers->end(); ++i) {
        if (!*i) continue;
        delete * i;
        *i = 0;
    }
}
template <class T> bool FinalTreeFiller::loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle, bool silent) {
    iEvent.getByLabel(tag, handle);
    if (!handle.isValid()) {
        if(!silent) edm::LogError("FinalTreeFiller") << "No valid handle found for '" << tag << "'!";
        return false;
    }
    return true;
}
