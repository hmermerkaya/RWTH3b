#ifndef ACHarvester_h
#define ACHarvester_h
/**
 Basic implementation of a harvesting step with the following abilities:
 \li print the contained instance of ACDataset if only one file is provided as commandline argument
 \li scan a single file for histograms and call a user defined function to process them
 \li combine histograms of one dataset from several (grid) jobs
 \li compare histograms of several dataset (if available) in a 'burst' mode (compare plot by plot). the user may overload this function to apply his prefered comparisons.
 \li do the final plots
 \li manipulate stored datasets

 Distinguish datasets by placing an instance of ACDataset in the top directory of the root files.
 \li split by dataset name in case all jobTypes are equal (comparison mode)
 \li split by jobType in case all dataset names are equal (validation mode)
 \li ignore datasets and merge all (merge mode)

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "boost/filesystem.hpp"

#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TChain.h"
#include "TTree.h"
#include <TKey.h>
#include <THashList.h>

#include "../../../ACDataFormats/ACGeneral/interface/ACDataset.h"

class ACHarvester {
public :
	ACHarvester(const std::string & outpath = ".", bool modifyInputFiles = false);
	virtual ~ACHarvester();
    
    /// recursivly read root files from directories and append them to fileNames
    void scanDirectories(const std::vector<std::string> & directories, std::vector<std::string> & fileNames);
    /// load root files and merge files from common datasets (separate datasets by analyzing the instance of ACDataset in the top directory of the root files, return success). set forcedMergeAll to merge all files ignoring their contained datasets
    bool loadRootFiles(const std::vector<std::string> & fileNames, bool printDatasets = false, bool forcedMergeAll = false);
    /// burst mode: compare histograms from different samples one by one
    void burstCompare();
    /// scan a single file for histograms and call user defined function processHistogram() to process them
    void scanFileForHistograms(const TFile * file);

protected:
    /// function called for every histogram found in the first dataset (to be defined by user, e.g. store it, ...)
    virtual void processHistogram(TH1 * h1, const std::string & path) {}
    /// merge histograms of several files (per default h1->Add(h2) is called), overload for more complex histograms
    virtual void mergeHistograms(TH1 * hist1, TH1 * hist2);
    /// compare histograms of several dataset one-by-one (can be replaced by user-defined function)
    virtual void compareHistograms(TH1 * hist1, TH1 * hist2, const ACDataset * dataset1, const ACDataset * dataset2, const std::string & storagePath, const std::string & path);
    /// try to recover mismatching binnings (e.g. different trigger menus), user-defined solution expected
    virtual bool recoverBinInconsistency(TH1 * h1, TH1 * h2) const;
    /// load user plot style for comparison plots
    virtual void compareStyle();
    /// persistantly modify dataset, only called if modifyInputFiles set to true. returns whether or not changes need to be stored
    virtual bool modifyDataset(ACDataset * dataset);

    /// whether or not to write changes to input datasets
    bool modifyInputFiles_;
    /// path to store the output
    std::string outpath_;
    /// list of merged files, one for each dataset
    std::vector<std::pair<ACDataset*, TFile*> > * mergedDatasets_;
    /// type of normalization applied in compareHistograms, 0: no norm, 1: sclae to common integral including outliers (default), 2: scale to lumi
    int normType_;
    
private:
    /// common label for all files (corresponds to either jobType or dataset name)
    std::string commonLabel_;
    /// operation mode: 1=compare mode, 2=validation mode, 3=merge mode
    int operationMode_;
    /// merge all files ignoring their contained datasets
    bool forcedMergeAll_;
    /// merged name stored into datasets when using forcedMergeAll_ (user defined during run time)
    std::string forcedMergeName_;
    
    /** read stored histograms and merge histograms according to mode, based on ROOT's hadd.C
     mode=0: combine files of one dataset from several (grid) jobs (arbitrary number of files)
     mode=1: compare histograms of different datasets (2 files)
     mode=2: process histograms of the first dataset
     */
    int mergeRootFiles(TDirectory * target, const TList * sourcelist, int count, int mode = 0);
    /** read stored objects and process each histogram according to user-definition in processHistogram().
        This function replaces the former mode==3 of mergeRootFiles (it caused a problem during deletion of the temporary created target)
     */
    void scanRootFile(const TDirectory * target, int depth);
    /// combine instances of ACDataset of equal id
    void combineDatasets(ACDataset * d1, const ACDataset * d2);
    
    /// recursivly read root files from directory and append them to fileNames
    void scanDirectory(const std::string & directory, std::vector<std::string> & fileNames, int level = 0);
    /// create required folders. storagePath is the path to the root file, path is the folder within
    bool createPlotDir(const std::string & storagePath, const std::string & path = "") const;
    /// prepare a writable dataset and call modifyDataset()
    bool prepareModificationOfDataset(TFile * file, TKey * key, ACDataset * dataset);
    /// test for equal binning
    bool checkCompatibleBinning(const TH1 * h1, const TH1 * h2 ) const;
    /// check that the axis limits of the histograms are the same
    bool checkAxisLimits(const TAxis *a1, const TAxis *a2 ) const;
    ///	create legend horizontally aligned on top of the pad
    void addHighLegend(TVirtualPad *can, const std::vector<TH1*> & histVect, const std::vector<std::string> & names);
    /// add a branding in comparison plots
    void addBranding(const std::string & branding);
    /// split a string into tokens separated by a specific delimiter (lazy copy from ACCommon.h)
    void tokenizePath(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ") const;
    /// helper function to read from command line input (std::cin)
    const std::string readInput() const;
};

#endif
