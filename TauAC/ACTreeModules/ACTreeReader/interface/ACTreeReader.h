#ifndef ACTreeReader_h
#define ACTreeReader_h
/**
 Class to access the content of the ACDataFormats.
 It requires a generic analyzer and provides the event loop.

 Based initially on Root's Make_Class.h

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <sys/stat.h>

#include "../../ACEvent/interface/ACEvent.h"
#include "../../ACProvenance/interface/ACProvenance.h"
#include "../../ACAnalyzer/interface/ACAnalyzer.h"

class ACTreeReader {
public:
    ACTreeReader() {};
    /// default constructor loading the entire list of available branches
    ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, TChain * chain = 0);
    /// default constructor loading the entire list of available branches (same as above) AND preparing a provenance object merged from another tree
    ACTreeReader(const std::vector<std::string>& fileNames, const std::string& treeName, const std::string& provenanceTreeName, TChain * chain = 0);
//    /// default constructor loading a user-defined list of branches (to be implemented later)
//    ACTreeReader(std::string fileName, std::string treeName, std::vector<std::string> branchNames);
    virtual ~ACTreeReader() {
        if (!fChain_) return;
        delete fChain_->GetCurrentFile();
    };

    /// performs the event loop and calls the analyze functions of the connected analyzer
    void loop(ACAnalyzer & analyzer, int maxEvents = -1); // WARNING!!! analyzer not const!!!
    /// number of entries in the chain
    int nentries() { return fChain_->GetEntriesFast(); }
    
    inline static bool GetVerbosity() { return verbosity_; }
    inline static void SetVerbosity(bool verbosity) { verbosity_ = verbosity; }
    
protected:
    /// list of filenames expected to contain the tree
    std::vector<std::string> fileNames_;
    /// list of branch names to be read from file
    std::vector<std::string> branchNames_;

    /// final event, which is read from the branches of the file. It will be provided within the event loop.
    ACEvent event_;
    /// final provenance, which is read from the branches of the file. It will be provided within the event loop.
    ACProvenance provenance_;

    ///event chain
    TChain *fChain_; //! pointer to the analyzed TChain
    ///provenance chain
    TChain *pChain_; //! pointer to the analyzed TChain
    Int_t fCurrent_; //! current Tree number in a TChain

    ///static verbosity switch
    static bool verbosity_;
    
    /// prepare files
    void loadFiles(const std::vector<std::string> & fileNames);
    /// connect all files to the chain
    TChain * setupChain(const std::string& treeName);
    /// connect all branches
    void Init(TChain *chain);
    /// connect all branches of the provenance
    void InitProvenance(TChain *chain);
    /// similar to loop, but here the elements of the chain are merged
    void mergeProvenance();
    /// Set the environment to read one entry
    Long64_t LoadTree(TChain * chain, Long64_t entry);
    Bool_t Notify();
    /// test the validity of the branch
    bool testBranch(TChain * chain, std::string branchname, int status, bool listAvailableBranches = false) const;
    /// print storage information
    void printTreeStatistics();
    /// print file input/output
    void printFileIO();
    /// split a string into tokens separated by a specific delimiter
    void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");    
};

#endif
