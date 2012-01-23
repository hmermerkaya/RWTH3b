#include "ACFrameworkModules/MCPileUpDistMaker/interface/MCPileUpDistMaker.h"

MCPileUpDistMaker::MCPileUpDistMaker(const edm::ParameterSet& iConfig):
pileupDistributions_(iConfig.getParameter<std::vector<edm::ParameterSet> >("pileUpDistributions")),
outputFileName_(iConfig.getParameter<std::string >("outputFile"))
{

}
MCPileUpDistMaker::~MCPileUpDistMaker() {}

void MCPileUpDistMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

}
void MCPileUpDistMaker::beginJob() {    
    TH1D * thisHist = 0;
    
    for (std::vector<edm::ParameterSet>::const_iterator idist = pileupDistributions_.begin(); idist != pileupDistributions_.end(); ++idist) {
        std::string thisName = idist->getParameter<std::string>("name");
        std::vector<double> thisDist = idist->getParameter<std::vector<double> >("dist");
        
        thisHist = fileService_->make<TH1D>(thisName.c_str(), thisName.c_str(), thisDist.size(), 0, thisDist.size()-1);
        unsigned int ibin = 1;
        for (std::vector<double>::const_iterator iter = thisDist.begin(); iter != thisDist.end(); ++iter, ibin++) {
            thisHist->SetBinContent(ibin, *iter);
        }
    }
}
void MCPileUpDistMaker::endJob() {

}
void MCPileUpDistMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

}
void MCPileUpDistMaker::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

}
void MCPileUpDistMaker::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {

}
void MCPileUpDistMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {

}

//define this as a plug-in
DEFINE_FWK_MODULE(MCPileUpDistMaker);
