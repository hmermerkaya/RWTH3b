#include "../interface/ACTreeReader.h"

ACTreeReader::ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, TChain *chain) {
    branchNames_.clear();
    loadFiles(fileNames);
    if (chain == 0) fChain_ = setupChain(treeName);
    Init(fChain_);
}
//ACTreeReader::ACTreeReader(std::string fileName, std::string treeName, std::vector<std::string> branchNames):
//ACTreeReader(fileName, treeName),
//branchNames_(branchNames)
//{
//}
ACTreeReader::ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, const std::string& provenanceTreeName, TChain *chain) {
    branchNames_.clear();
    loadFiles(fileNames);
    if (chain == 0) fChain_ = setupChain(treeName);
    Init(fChain_);

    pChain_ = setupChain(provenanceTreeName);
    InitProvenance(pChain_);
    mergeProvenance();
}

void ACTreeReader::loadFiles(const std::vector<std::string> & fileNames) {
    //check for white spaces and split the string if needed
    for (std::vector<std::string>::const_iterator ifile = fileNames.begin(); ifile != fileNames.end(); ++ifile) {
        std::vector<std::string> tmp;
        tokenize(*ifile, tmp);
        fileNames_.insert(fileNames_.end(), tmp.begin(), tmp.end());
    }
}
TChain * ACTreeReader::setupChain(const std::string& treeName) {
    TChain * chain = new TChain(treeName.c_str());
    
    for (std::vector<std::string>::const_iterator ifile = fileNames_.begin(); ifile != fileNames_.end(); ++ifile) {
        std::string abspath = "";
        if (ifile->substr(0,1) != "/" && ifile->substr(0,4) != "dcap") {
            abspath = "file:"+std::string(gSystem->pwd())+"/"+(*ifile);
        } else if (ifile->substr(0,4) == "dcap") {
            abspath = (*ifile);
        } else {
            abspath = "file:"+(*ifile);
        }
        int status = chain->AddFile(abspath.c_str(),0);
        if (status != 1) {
            std::cerr << "ERROR: File '" << abspath << "' could not be connected to the chain of input files! Maybe the tree '" << treeName << "' does not exist in this file?!" << std::endl;
            throw 404;
        }
    }
    return chain;
}

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
    //fChain_->AddBranchToCache("*",kTRUE);//speedup (gain like 1s in 10000 events)
    
    Long64_t nbytes = 0, nb = 0;
    int run = -1;
    int lumi = -1;
    fCurrent_ = -1;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(fChain_, jentry);
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
            analyzer.beginJob(provenance_);//empty object if there is no pChain_
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
void ACTreeReader::mergeProvenance() {
    if (pChain_ == 0) return;
    Long64_t nentries = pChain_->GetEntriesFast();
    
    //increase cache size
    int cachesize = 30000000;   //30 MBytes
    pChain_->SetCacheSize(cachesize);
    pChain_->SetCacheLearnEntries();
    //pChain_->AddBranchToCache("*",kTRUE);//speedup (gain like 1s in 10000 events)
    
    fCurrent_ = -1;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(pChain_, jentry);
        if (ientry < 0) break;
        pChain_->GetEntry(jentry);
        // now merge the provenance contents
        provenance_.merge();
    }
}

void ACTreeReader::Init(TChain *chain) {
    /// Set branch addresses and branch pointers
    if (!chain) return;
    
    //    chain->SetBranchStatus("*",0);  // disable all branches
    //    chain->SetBranchStatus("branchname",1);  // activate branchname
    event_ = ACEvent(); // empty default c'tor. fill by set functions
    
    int status = 0;
    
    bool listAvailableBranches = true;//report only once
    
    status = chain->SetBranchAddress("ACEventInfo", event_.linkEventInfo());
    listAvailableBranches = testBranch(chain, "ACEventInfo", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACEventWeight", event_.linkEventWeight());
    listAvailableBranches = testBranch(chain, "ACEventWeight", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACEventGlobals", event_.linkEventGlobals());
    listAvailableBranches = testBranch(chain, "ACEventGlobals", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACTrigger", event_.linkTrigger());
    listAvailableBranches = testBranch(chain, "ACTrigger", status, listAvailableBranches);

    status = chain->SetBranchAddress("triggerObjects", event_.linkTriggerObjects());
    listAvailableBranches = testBranch(chain, "triggerObjects", status, listAvailableBranches);

    status = chain->SetBranchAddress("ACBeamSpot", event_.linkBeamSpot());
    listAvailableBranches = testBranch(chain, "ACBeamSpot", status, listAvailableBranches);

    status = chain->SetBranchAddress("offlinePV", event_.linkOfflinePV());
    listAvailableBranches = testBranch(chain, "offlinePV", status, listAvailableBranches);
    status = chain->SetBranchAddress("reducedPV", event_.linkReducedPV());
    listAvailableBranches = testBranch(chain, "reducedPV", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("generator", event_.linkGenerator());
    listAvailableBranches = testBranch(chain, "generator", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("genTauDecays", event_.linkGenTauDecays());
    listAvailableBranches = testBranch(chain, "genTauDecays", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("Tracks", event_.linkTracks());
    listAvailableBranches = testBranch(chain, "Tracks", status, listAvailableBranches);

    status = chain->SetBranchAddress("Muons", event_.linkMuons());
    listAvailableBranches = testBranch(chain, "Muons", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("Electrons", event_.linkElectrons());
    listAvailableBranches = testBranch(chain, "Electrons", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACFittedThreeProngParticles", event_.linkFittedThreeProngParticles());
    listAvailableBranches = testBranch(chain, "ACFittedThreeProngParticles", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACFittedThreeProngDecays", event_.linkTauDecays());
    listAvailableBranches = testBranch(chain, "ACFittedThreeProngDecays", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("PFJets", event_.linkPFJets());
    listAvailableBranches = testBranch(chain, "PFJets", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("PFTaus", event_.linkPFTaus());
    listAvailableBranches = testBranch(chain, "PFTaus", status, listAvailableBranches);
    
    status = chain->SetBranchAddress("ACPileupInfo", event_.linkPileup());
    listAvailableBranches = testBranch(chain, "ACPileupInfo", status, listAvailableBranches);
}
void ACTreeReader::InitProvenance(TChain *chain) {
    /// Set branch addresses and branch pointers
    if (!chain) return;
    
    int status = 0;
    
    status = chain->SetBranchAddress("ModuleCounter", provenance_.linkCounter());
    testBranch(chain, "ModuleCounter", status, false);
    
    status = chain->SetBranchAddress("ModuleCounterWeighted", provenance_.linkCounterWeighted());
    testBranch(chain, "ModuleCounterWeighted", status, false);
    
}

Long64_t ACTreeReader::LoadTree(TChain * chain, Long64_t entry) {
    /// Set the environment to read one entry
    if (!chain) return -5;
    Long64_t centry = chain->LoadTree(entry);
    if (centry < 0) return centry;
    if (!chain->InheritsFrom(TChain::Class()))  return centry;
    TChain *newchain = (TChain*)chain;
    if (newchain->GetTreeNumber() != fCurrent_) {
        fCurrent_ = newchain->GetTreeNumber();
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
bool ACTreeReader::testBranch(TChain * chain, std::string branchname, int status, bool listAvailableBranches) const {
    //std::cout<<"branchstatus "<<branchname<<": "<<status<<std::endl;
    if (!chain->GetBranchStatus(branchname.c_str())) {
        std::cout<<"ACTreeReader::checkBranch: branch '"<<branchname<<"' is invalid and will be ignored!"<<std::endl;
        //chain->SetBranchStatus(branchname.c_str(), 0);
        if (listAvailableBranches) {
            std::cout<<"Here is a list of available branches:"<<std::endl;
            chain->GetListOfBranches()->Print();
        }
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
