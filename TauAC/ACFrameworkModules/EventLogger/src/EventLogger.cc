#include "../interface/EventLogger.h"

EventLogger::EventLogger(const edm::ParameterSet& iConfig):
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
pathName_(iConfig.getParameter<std::string>("pathName")),
modules_(iConfig.getUntrackedParameter<std::vector<std::string> >("modules")),
pileUpDistributionFileMC_(iConfig.getUntrackedParameter<std::string>("pileUpDistributionFileMC", "")),
pileUpDistributionHistMC_(iConfig.getUntrackedParameter<std::string>("pileUpDistributionHistMC", "")),
pileUpDistributionFileData_(iConfig.getUntrackedParameter<std::string>("pileUpDistributionFileData", "")),
pileUpDistributionHistData_(iConfig.getUntrackedParameter<std::string>("pileUpDistributionHistData", ""))
{
    jobTree_ = 0;
    lumiWeights_ = 0;
}

EventLogger::~EventLogger() {
}

// ------------ method called on each new Event  ------------
void EventLogger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    //store the pass or fail result of each path (reads TriggerResults of the current process (not HLT))
    
    // calculate event weights in case of MC
    double eventWeight = 1.0;
    if (lumiWeights_) {
        const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
        eventWeight = lumiWeights_->weight3D(*iEventB);
    }

    std::string allEvtLabel = "event_0_AllEvents_cntSuccess";
    // first event setup
    if (labels_.size() == 0) {
        //make sure to be first element in map. stick to scheme 'event_orderID_name_type'
        labels_.push_back(allEvtLabel);
        filterResults_.add(allEvtLabel);
        if (lumiWeights_) filterResultsWeighted_.add(allEvtLabel);
    }
    filterResults_.count(allEvtLabel);
    if (lumiWeights_) filterResultsWeighted_.count(allEvtLabel, eventWeight);

    if (!HLTCP_.inited()) {
        std::cout<<"EventLogger::filter: HLT config provider not initiated!"<<std::endl;
        return;
    }
    
    edm::Handle<edm::TriggerResults> pathFilterResults;
    if (loadCollection(iEvent, triggerResultsTag_, pathFilterResults)) {
        edm::TriggerNames triggerNames = iEvent.triggerNames(*pathFilterResults);
        std::map<std::string, std::vector<std::string> > moduleNamesPerPath;
        for (unsigned int i = 0; i != triggerNames.triggerNames().size(); ++i) {
            const std::string & pathName = triggerNames.triggerName(i);
            if (pathName != pathName_) continue;
            
            std::vector<std::string> modules = HLTCP_.moduleLabels(pathName);//label(s) of module(s) in a trigger path
            moduleNamesPerPath[pathName] = modules;
//            std::cout<<"\twork on path "<<pathName<<"("<<pathFilterResults->accept(i)<<") with "<<modules.size()<<" filter modules."<<std::endl;
//            std::cout<<"\tIndex of the failed module is "<<pathFilterResults->index(i)<<std::endl;
//            std::cout<<"\tName of the failed module is "<<modules.at(pathFilterResults->index(i))<<std::endl;
            
            //prepare the storage once
            if (filterResults_.counts().size() == 1) {
                unsigned int index = 0;
                for(std::vector<std::string>::const_iterator module = modules.begin(); module != modules.end(); ++module, index++) {
                    if (skipModule(*module)) continue;
                    //orderID is the module index starting from 1(!) to create slot for all event counter.
                    //stick to scheme 'event_orderID_name_type'
                    std::string pattern = "event_" + any2string(index+1) + "_" + *module + "_cntSuccess";
                    filterResults_.add(pattern);
                    if (lumiWeights_) filterResultsWeighted_.add(pattern);
                    labels_.push_back(pattern);
                }
            }
            unsigned int index = 0;
            for(std::vector<std::string>::const_iterator module = modules.begin(); module != modules.end(); ++module, index++) {
                if (skipModule(*module)) continue;
                //orderID is the module index starting from 1(!) to create slot for all event counter.
                //stick to scheme 'event_orderID_name_type'
                std::string pattern = "event_" + any2string(index+1) + "_" + *module + "_cntSuccess";
                if (pathFilterResults->accept(i)) {
                    filterResults_.count(pattern);
                    if (lumiWeights_) filterResultsWeighted_.count(pattern, eventWeight);
                } else {
                    //only count successful modules (if accept the last module would define the decision)
                    if (index < pathFilterResults->index(i)) {
                        filterResults_.count(pattern);
                        if (lumiWeights_) filterResultsWeighted_.count(pattern, eventWeight);
                    }
                }
            }
        }
    }
}

void EventLogger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // initialize HLTConfigProvider
    bool changed = false;//paths may not change
    HLTCP_.init(iRun, iSetup, triggerResultsTag_.process(), changed);
}

// ------------ method called once each job just before starting event loop  ------------
void EventLogger::beginJob() {
    //edm::LogVerbatim("EventLogger")<<"--> [EventLogger] is logging";
    jobTree_ = fileService_->make<TTree>("ModuleSummary", (std::string("Summary of all modules in path ") + pathName_).c_str());
    jobTree_->Branch("ModuleCounter", &filterResults_, 32000, 0);
    
    // initialize pileup reweighting (it might be useful to renew the initialization for every run!?)
    if (pileUpDistributionFileMC_ != "" && pileUpDistributionFileData_ != "") {
        lumiWeights_ = new edm::Lumi3DReWeighting(pileUpDistributionFileMC_, pileUpDistributionFileData_, pileUpDistributionHistMC_, pileUpDistributionHistData_, "");
        lumiWeights_->weight3D_init(1.0);
        jobTree_->Branch("ModuleCounterWeighted", &filterResultsWeighted_, 32000, 0);
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void EventLogger::endJob() {
    jobTree_->Fill();
    edm::LogVerbatim("EventLogger")<<"--> [EventLogger]";//<<"\n"<<filterResults_;
    edm::LogVerbatim("EventLogger")<<printCounter(filterResults_);
    if (lumiWeights_) edm::LogVerbatim("EventLogger")<<"--> [EventLoggerWeighted]\n"<<printCounter(filterResultsWeighted_);
    
    delete lumiWeights_;
}

const std::string EventLogger::stripCommonPattern(const std::string & label) const {
    // expect label = event_orderID_moduleName_cntSuccess
    if (label.size() < 18) return label;
    if (label.substr(0, 6) != "event_") return label;
    if (label.substr(label.size() - 11) != "_cntSuccess") return label;

    return label.substr(6, label.size() - 11 - 6);
}
const std::string EventLogger::printCounter(const ACCounter & filterResults) const {
    std::ostringstream buffer;
    
    //need the order here
    unsigned int labelWidth = 0;
    for (std::vector<std::string>::const_iterator label = labels_.begin(); label != labels_.end(); ++label) {
        std::string strippedLabel = "[" + stripCommonPattern(*label) + "]";
        labelWidth = std::max(labelWidth, (unsigned int)strippedLabel.size());
    }
    unsigned int evtWidth = 0;
    std::vector<std::string>::const_iterator lastLabel = labels_.begin();
    for (std::vector<std::string>::const_iterator label = labels_.begin(); label != labels_.end(); ++label) {
        std::map<std::string, double>::const_iterator module = filterResults.counts().find(*label);
        std::map<std::string, double>::const_iterator lastModule = filterResults.counts().find(*lastLabel);
        //        std::map<std::string, double>::const_iterator moduleW = filterResultsWeighted_.counts().find(*label);
        //        std::map<std::string, double>::const_iterator lastModuleW = filterResultsWeighted_.counts().find(*lastLabel);
        std::string strippedLabel = "[" + stripCommonPattern(module->first) + "]";
        if (label != lastLabel) {
            ++lastLabel;
            std::string postfix = "";
            if (strippedLabel.find("_-")!=std::string::npos) postfix = "\t(filter ignored)";
            buffer<<"    "<<std::setw(labelWidth)<<strippedLabel<<" Efficiency: "<<std::setw(evtWidth)<<module->second<<"/"<<std::setw(evtWidth)<<lastModule->second<<" = "<<std::setw(4)<<filterResults.ratio(module->first, lastModule->first, -3) * 100<<" %"<<postfix<<std::endl;
            //<<" weighted: "<<
        } else {
            evtWidth = any2string(module->second).size();
            //all event counter
            buffer<<"    "<<std::setw(labelWidth)<<strippedLabel<<" Counter:    "<<std::setw(2*evtWidth+1)<<module->second<<std::endl;
        }
    }
    if (labels_.size() > 2) {
        buffer<<"    "<<std::string(labelWidth + std::string(" Efficiency: ").size() + 2*evtWidth+1 + 9,'-')<<std::endl;
        std::map<std::string, double>::const_iterator firstModule = filterResults.counts().find(labels_.front());
        std::map<std::string, double>::const_iterator lastModule = filterResults.counts().find(labels_.back());
        buffer<<"    "<<std::setw(labelWidth)<<"[Summary]"<<" Efficiency: "<<std::setw(evtWidth)<<lastModule->second<<"/"<<std::setw(evtWidth)<<firstModule->second<<" = "<<std::setw(4)<<filterResults.ratio(lastModule->first, firstModule->first, -3) * 100<<" %"<<std::endl;
    }

    return buffer.str();
}

bool EventLogger::skipModule(const std::string & module) const {
    // do not skip any module if the module list is empty
    if (modules_.size() == 0) return false;
    // skip possible empty entries (never seen yet)
    if (module.size() < 1) return true;
    
    std::string moduleStripped = module;
    if (module.substr(0,1)=="-") moduleStripped = module.substr(1);//deactivated EDFilter carry a leading minus in their name
    for (std::vector<std::string>::const_iterator label = modules_.begin(); label != modules_.end(); ++label) {
        // do not skip the module if it's in the module list
        if (*label == moduleStripped) return false;
    }
    
    return true;
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventLogger);
