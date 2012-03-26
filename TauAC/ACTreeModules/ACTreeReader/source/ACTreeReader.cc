#include "../interface/ACTreeReader.h"

ACTreeReader::ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, TChain *chain) {
    branchNames_.clear();
    
    if (chain == 0) {
        fChain_ = new TChain(treeName.c_str());
        
        //check for white spaces and split the string if needed
        for (std::vector<std::string>::const_iterator ifile = fileNames.begin(); ifile != fileNames.end(); ++ifile) {
            std::vector<std::string> tmp;
            tokenize(*ifile, tmp);
            fileNames_.insert(fileNames_.end(), tmp.begin(), tmp.end());
        }
        for (std::vector<std::string>::const_iterator ifile = fileNames_.begin(); ifile != fileNames_.end(); ++ifile) {
            std::string abspath = "";
            if (ifile->substr(0,1) != "/" && ifile->substr(0,4) != "dcap") {
                abspath = "file:"+std::string(gSystem->pwd())+"/"+(*ifile);
            } else if (ifile->substr(0,4) == "dcap") {
                abspath = (*ifile);
            } else {
                abspath = "file:"+(*ifile);
            }
            int status = fChain_->AddFile(abspath.c_str(),0);
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

bool ACTreeReader::verbosity_ = false; 

void ACTreeReader::loop(ACAnalyzer & analyzer, int maxEvents) {
    if (fChain_ == 0) return;
    Long64_t nentries = fChain_->GetEntriesFast();
    if (maxEvents!=-1 && maxEvents<nentries) {
        if (verbosity_) printf("ACTreeReader::Loop: analyze only %d events of %d events.\n", (int)maxEvents, (int)nentries);
        nentries = maxEvents;
    } else {
        if (verbosity_) printf("ACTreeReader::Loop: analyze %d events.\n", (int)nentries);
    }
    
    //increase cache size
    int cachesize = 30000000;   //30 MBytes
    fChain_->SetCacheSize(cachesize);
    fChain_->SetCacheLearnEntries();
    
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
            if (verbosity_) printTreeStatistics();
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
            if (verbosity_) printFileIO();
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
    
    status = fChain_->SetBranchAddress("ACEventWeight", event_.linkEventWeight());
    testBranch("ACEventWeight", status);
    
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
    
    status = fChain_->SetBranchAddress("Tracks", event_.linkTracks());
    testBranch("Tracks", status);

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
void ACTreeReader::printTreeStatistics() {
    float totBytes = fChain_->GetTree()->GetTotBytes();
    float zipBytes = fChain_->GetTree()->GetZipBytes();
    printf("ACTreeReader::printTreeStatistics:\n total size of tree '%s' (raw/compressed) in MB: %.1f / %.1f\n", fChain_->GetTree()->GetName(), totBytes/1024.0/1024.0, zipBytes/1024.0/1024.0);
    unsigned int numberOfBranches = fChain_->GetTree()->GetListOfBranches()->GetEntriesFast();
    
    for (unsigned int i = 0; i < numberOfBranches; i++) {
        TBranch * thisBranch = (TBranch*)(fChain_->GetTree()->GetListOfBranches()->At(i));
        float thisTotBytes = thisBranch->GetTotBytes();
        float thisZipBytes = thisBranch->GetZipBytes();
        printf("  branch '%s': size (raw/compressed) in MB: %.1f / %.1f --> %.1f%% / %.1f%%\n", thisBranch->GetName(), thisTotBytes/1024.0/1024.0, thisZipBytes/1024.0/1024.0, thisTotBytes/totBytes*100.0, thisZipBytes/zipBytes*100.0);
    }
}
void ACTreeReader::printFileIO() {
    double bytesRead = 0.0;
    double readCalls = 0.0;
    //    TObjArray * listOfFiles = fChain_->GetListOfFiles();
    //    unsigned int numberOfFiles = listOfFiles->GetEntriesFast();
    //    for (unsigned int i = 0; i < numberOfFiles; i++) {
    //        TFile * thisFile = (TFile*)(listOfFiles->At(i));
    //        bytesRead += thisFile->GetBytesRead();
    //        readCalls += thisFile->GetReadCalls();
    //    }
    bytesRead = fChain_->GetCurrentFile()->GetBytesRead();
    readCalls = fChain_->GetCurrentFile()->GetReadCalls();
    printf("ACTreeReader::printFileIO: reading %.1f MB in %d transactions --> %.3f MB/call\n",bytesRead/1024.0/1024.0, int(readCalls), bytesRead/1024.0/1024.0/readCalls);
}
void ACTreeReader::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
    /// Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    /// Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos) {
        /// Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        /// Skip delimiters. Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        /// Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
