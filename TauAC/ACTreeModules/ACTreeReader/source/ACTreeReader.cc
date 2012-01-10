#include "../interface/ACTreeReader.h"

ACTreeReader::ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, TChain *chain):
fileNames_(fileNames)
{
    branchNames_.clear();

    if (chain == 0) {
        fChain_ = new TChain(treeName.c_str());
        for (std::vector<std::string>::const_iterator ifile = fileNames_.begin(); ifile != fileNames_.end(); ++ifile) {
            std::string abspath = ""; 
            if (ifile->substr(0,1) != "/") abspath = std::string(gSystem->pwd())+"/";
            int status = fChain_->AddFile(("file:"+abspath+(*ifile)).c_str(),0);
            if (status != 1) {
                std::cerr << "ERROR: File '" << abspath << "' could not be connected to the chain of input files! Maybe the tree '" << treeName << "' does not exist in this file?!" << std::endl;
            }
        }
    }
    Init(fChain_);
}

//ACTreeReader::ACTreeReader(std::string fileName, std::string treeName, std::vector<std::string> branchNames):
//ACTreeReader(fileName, treeName),
//branchNames_(branchNames)
//{
//}

void ACTreeReader::loop(ACAnalyzer & analyzer, int maxEvents) {
    if (fChain_ == 0) return;
    Long64_t nentries = fChain_->GetEntriesFast();
    if (maxEvents!=-1 && maxEvents<nentries) {
        printf("ACTreeReader::Loop: analyze only %d events of %d events.\n", (int)maxEvents, (int)nentries);
        nentries = maxEvents;
    } else {
        printf("ACTreeReader::Loop: analyze %d events.\n", (int)nentries);
    }

    Long64_t nbytes = 0, nb = 0;
    int run = -1;
    int lumi = -1;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain_->GetEntry(jentry);
        nbytes += nb;
//        if (nentries < 11) printf("evt %d:\n", (int)jentry+1);
//        if (nentries > 10 && nentries < 101)
//            if ((jentry+1) % (int)round(nentries/10) == 0 || jentry==nentries-1) {
//                printf("evt %d:\r", (int)jentry+1);
//                fflush(stdout);
//            }
//        if (nentries > 100)
//            if ((jentry+1) % (int)round(nentries/100) == 0 || jentry==nentries-1) {
//                printf("evt %d:\r", (int)jentry+1);
//                fflush(stdout);
//            }
//        if (jentry==nentries-1) printf("\n");

        int irun = event_.eventInfo()->runID();
        int ilumi = event_.eventInfo()->lumiID();

        if (jentry == 0) {
            analyzer.beginJob();
            run = irun;
            analyzer.beginRun();
            lumi = ilumi;
            analyzer.beginLuminosityBlock();
        }

        if (lumi != ilumi) {
            analyzer.endLuminosityBlock();
            lumi = ilumi;
            analyzer.beginLuminosityBlock();
        }

        if (run != irun) {
            analyzer.beginRun();
            run = irun;
            analyzer.endRun();
        }

        analyzer.analyze(event_);

        if (jentry == nentries-1) {
            analyzer.endLuminosityBlock();
            analyzer.endRun();
            analyzer.endJob();
        }
    }
}

void ACTreeReader::Init(TChain *chain) {
    /// Set branch addresses and branch pointers
    if (!chain) return;
    fCurrent_ = -1;


//    fChain_->SetBranchStatus("*",0);  // disable all branches
//    fChain_->SetBranchStatus("branchname",1);  // activate branchname
    event_ = ACEvent(); // empty default c'tor. fill by set functions

    int status = 0;

    status = fChain_->SetBranchAddress("ACEventInfo", event_.linkEventInfo());
    testBranch("ACEventInfo", status);

    status = fChain_->SetBranchAddress("ACEventGlobals", event_.linkEventGlobals());
    testBranch("ACEventGlobals", status);

    status = fChain_->SetBranchAddress("ACTrigger", event_.linkTrigger());
    testBranch("ACTrigger", status);

    status = fChain_->SetBranchAddress("offlinePV", event_.linkOfflinePV());
    testBranch("offlinePV", status);
    status = fChain_->SetBranchAddress("reducedPV", event_.linkReducedPV());
    testBranch("reducedPV", status);

    status = fChain_->SetBranchAddress("generator", event_.linkGenerator());
    testBranch("generator", status);

    status = fChain_->SetBranchAddress("genTauDecays", event_.linkGenTauDecays());
    testBranch("genTauDecays", status);

    status = fChain_->SetBranchAddress("Muons", event_.linkMuons());
    testBranch("Muons", status);

    status = fChain_->SetBranchAddress("Electrons", event_.linkElectrons());
    testBranch("Electrons", status);

    status = fChain_->SetBranchAddress("ACFittedThreeProngParticles", event_.linkFittedThreeProngParticles());
    testBranch("ACFittedThreeProngParticles", status);

    status = fChain_->SetBranchAddress("ACFittedThreeProngDecays", event_.linkTauDecays());
    testBranch("ACFittedThreeProngDecays", status);

    status = fChain_->SetBranchAddress("PFJets", event_.linkPFJets());
    testBranch("PFJets", status);

    status = fChain_->SetBranchAddress("PFTaus", event_.linkPFTaus());
    testBranch("PFTaus", status);

    status = fChain_->SetBranchAddress("ACPileupInfo", event_.linkPileup());
    testBranch("ACPileupInfo", status);
}

Long64_t ACTreeReader::LoadTree(Long64_t entry) {
    /// Set the environment to read one entry
    if (!fChain_) return -5;
    Long64_t centry = fChain_->LoadTree(entry);
    if (centry < 0) return centry;
    if (!fChain_->InheritsFrom(TChain::Class()))  return centry;
    TChain *chain = (TChain*)fChain_;
    if (chain->GetTreeNumber() != fCurrent_) {
        fCurrent_ = chain->GetTreeNumber();
        Notify();
    }
    return centry;
}
Bool_t ACTreeReader::Notify() {
    /// The Notify() function is called when a new file is opened. This
    /// can be either for a new TTree in a TChain or when when a new TTree
    /// is started when using PROOF. It is normally not necessary to make changes
    /// to the generated code, but the routine can be extended by the
    /// user if needed. The return value is currently not used.

    return kTRUE;
}
bool ACTreeReader::testBranch(std::string branchname, int status) const {
    if (!fChain_->GetBranchStatus(branchname.c_str())) {
        printf("ACTreeReader::checkBranch: branch %s is invalid!\n", branchname.c_str());
        return false;
    }
    return true;
}
