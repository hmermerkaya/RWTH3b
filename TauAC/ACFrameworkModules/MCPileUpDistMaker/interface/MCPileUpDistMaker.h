/**
 CMSSW analyzer to store a root file with MC pileup histograms.

 @author Lars Perchalla & Philip Sauerland
 @date 2012
 */
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class MCPileUpDistMaker : public edm::EDAnalyzer {
public:
    explicit MCPileUpDistMaker(const edm::ParameterSet&);
    ~MCPileUpDistMaker();

private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    /// framework file service
    edm::Service<TFileService> fileService_;
    /// input tags
    std::vector<edm::ParameterSet> pileupDistributions_;
    std::string outputFileName_;
};
