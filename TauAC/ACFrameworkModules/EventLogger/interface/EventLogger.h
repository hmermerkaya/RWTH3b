/**
 CMSSW analyzer logging status of all paths derived from the TriggerResults object stored for each path.
 The TriggerResultInserter is run after all normal paths have completed, but before any end paths have run.
 Therefore this module has to be called from an EndPath.

 @author Lars Perchalla
 @date 2013
 */
#include <memory>

#include <TTree.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ACDataFormats/ACGeneral/interface/ACCounter.h"

#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

class EventLogger : public edm::EDAnalyzer {
public:
    explicit EventLogger(const edm::ParameterSet&);
    ~EventLogger();

private:
    virtual void beginJob() ;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    /// save way to fill a collection from a inputtag into a given handle
    template <class T> bool loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle, bool silent = false);
    ///convert any to string
    template <class T> std::string any2string(T i) const;
    /// strip common part of pattern event_orderID_moduleName_cntSuccess from label
    const std::string stripCommonPattern(const std::string & label) const;
    /// print the counters
    const std::string printCounter(const ACCounter & filterResults) const;
    /// test the current module should be monitored or not
    bool skipModule(const std::string & module) const;
    
    /// framework file service
    edm::Service<TFileService> fileService_;
    /// filled only once, but need tree to store own objects with TFileService
    TTree * jobTree_;

    /// TriggerResults object to monitor
    edm::InputTag triggerResultsTag_;
    /// path to monitor
    std::string pathName_;
    /// module labels to monitor (if empty, log all modules)
    std::vector<std::string> modules_;
    /// pileup information for event weighting
    std::string pileUpDistributionFileMC_, pileUpDistributionHistMC_, pileUpDistributionFileData_, pileUpDistributionHistData_;
    HLTConfigProvider HLTCP_;
    edm::Lumi3DReWeighting * lumiWeights_;
    
    /// ordered list of counters (same size as filterResults_)
    std::vector<std::string> labels_;
    /// log the status of all modules (use scheme event_orderID_moduleName_cntSuccess) in the current path and count their success (paths may not change). orderID is the module index starting from 1(!) to create slot for all event counter.
    ACCounter filterResults_;
    /// same as filterResults_ but including event weights
    ACCounter filterResultsWeighted_;
};

template <class T> std::string EventLogger::any2string(T i) const {
    std::ostringstream buffer;
    buffer << i;
    return buffer.str();
}

template <class T> bool EventLogger::loadCollection(const edm::Event& iEvent, const edm::InputTag & tag, edm::Handle<T>& handle, bool silent) {
    iEvent.getByLabel(tag, handle);
    if (!handle.isValid()) {
        if(!silent) edm::LogError("EventLogger") << "No valid handle found for '" << tag << "'!";
        return false;
    }
    return true;
}
