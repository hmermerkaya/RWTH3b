#define ACHarvester_cxx
#include "../interface/ACHarvester.h"

ACHarvester::ACHarvester(const std::string & outpath, bool modifyInputFiles) {
    outpath_ = outpath;
    mergedDatasets_ = new std::vector<std::pair<ACDataset*, TFile*> >;
    normType_ = 1;
    commonLabel_ = "";
    operationMode_ = 0;
    modifyInputFiles_ = modifyInputFiles;
    forcedMergeAll_ = false;
    forcedMergeName_ = "";
}
ACHarvester::~ACHarvester() {
    for (std::vector<std::pair<ACDataset*, TFile*> >::iterator i = mergedDatasets_->begin(); i < mergedDatasets_->end(); ++i) {
        if (i->first) {
            delete i->first;
        }
        if (i->second) {
            //std::cout<<"cleanup "<<i->second->GetName()<<std::endl;
            i->second->Close();
            delete i->second;
        }
    }
    delete mergedDatasets_;
    //std::cout<<"d'tor ACHarvester"<<std::endl;
}

// main functions
void ACHarvester::scanDirectories(const std::vector<std::string> & directories, std::vector<std::string> & fileNames) {
    for (std::vector<std::string>::const_iterator dir = directories.begin(); dir != directories.end(); ++dir) {
        scanDirectory(*dir, fileNames);
        //std::cout<<"found "<<fileNames.size()<<std::endl;
    }
}
bool ACHarvester::loadRootFiles(const std::vector<std::string> & fileNames, bool printDatasets, bool forcedMergeAll) {
    if (fileNames.size() < 1) {
        std::cout<<"ACHarvester::loadRootFiles: ERROR! Empty file list!"<<std::endl;
        return false;
    }
    forcedMergeAll_ = forcedMergeAll;
    std::vector<std::pair<ACDataset*, TFile*> > allDatasets;//one for each input file
    std::string jobType = "";
    std::string name = "";
    bool splitByDatasetName = false;
    bool splitByJobType = false;
    bool modified = false;
    for (std::vector<std::string>::const_iterator ifile = fileNames.begin(); ifile != fileNames.end(); ++ifile) {
        std::string abspath = "";
        if (ifile->substr(0,1) != "/" && ifile->substr(0,4) != "dcap") {
            abspath = "file:"+std::string(gSystem->pwd())+"/"+(*ifile);
        } else if (ifile->substr(0,4) == "dcap") {
            abspath = (*ifile);
        } else {
            abspath = "file:"+(*ifile);
        }
        TFile * file = 0;
        if (!modifyInputFiles_) {
            file = TFile::Open(abspath.c_str(), "READ");
        } else {
            std::cout<<"ACHarvester::loadRootFiles: WARNING! Opening file in write mode at "<<*ifile<<std::endl;
            file = TFile::Open(abspath.c_str(), "UPDATE");
        }
        if (!file) {
            std::cout<<"ACHarvester::loadRootFiles: ERROR! No file at "<<*ifile<<std::endl;
            return false;
        }
        
        // scan for instance of ACDataset in top folder
        ACDataset * dataset = 0;
        file->cd();
        TDirectory *current_sourcedir = gDirectory;
        Bool_t status = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kFALSE);
        TIter nextkey(current_sourcedir->GetListOfKeys());
        TKey *key, *oldkey=0;
        while ((key = (TKey*)nextkey())) {
            //keep only the highest cycle number for each key
            if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
            TObject *obj = key->ReadObj();
            if (obj->IsA()->InheritsFrom("ACDataset")) {
                if (!dataset) {
                    dataset = (ACDataset*)obj;
                    modified = prepareModificationOfDataset(file, key, dataset);
                    if (printDatasets) {
                        std::cout<<*dataset<<std::endl;
                    }
                    if (name=="" && jobType=="") {
                        name = dataset->name();
                        jobType = dataset->jobType();
                    } else {
                        if (dataset->name() != name) {
                            //detected different names
                            splitByDatasetName = true;
                        }
                        if (dataset->jobType() != jobType) {
                            //detected different types
                            splitByJobType = true;
                        }
                    }
                    if (modified) {
                        break;
                    }
                } else {
                    std::cout<<"ACHarvester::identifyRootFiles:ERROR! multiple datasets found in "<<ifile->c_str()<<std::endl;
                    return 101;
                }
            } else {
                continue;
            }
        }
        TH1::AddDirectory(status);
        if (!dataset) {
            std::cout<<"ACHarvester::identifyRootFiles:ERROR! no dataset found in "<<ifile->c_str()<<std::endl;
            return 101;
        }
        allDatasets.push_back(std::make_pair(dataset, file));
    }
    if (modified) {
        std::cout<<"ACHarvester::loadRootFiles:Input files have been modified. Skip further harvesting."<<std::endl;
        return false;
    }
    if (!forcedMergeAll_) {
        if (splitByJobType && splitByDatasetName) {
            std::cout<<"ACHarvester::loadRootFiles:ERROR! Mismatch both in jobTypes of dataset names. Harvesting impossible."<<std::endl;
            throw 404;
        }
        if (splitByJobType && !splitByDatasetName) {
            operationMode_ = 2;
            commonLabel_ = name;
            std::cout<<"--> Validation mode: Compare datasets of different jobTypes and common name "<<commonLabel_<<"."<<std::endl;
        }
        if (!splitByJobType && splitByDatasetName) {
            operationMode_ = 1;
            commonLabel_ = jobType;
            std::cout<<"--> Comparison mode: Compare datasets of different names and common jobType "<<commonLabel_<<"."<<std::endl;
        }
        if (!splitByJobType && !splitByDatasetName) {
            operationMode_ = 1;
            commonLabel_ = jobType;
            if (fileNames.size()>1) std::cout<<"--> Burst mode: Datasets are equal in jobType and dataset name. Merge all of them."<<std::endl;
        }
    } else {
        operationMode_ = 3;
        // try to find common jobType
        if (!splitByJobType) {
            commonLabel_ = jobType;
            std::cout<<"--> Forced burst mode: Datasets names have been ignored. Merge all files of common jobType "<<jobType<<"."<<std::endl;
        } else {
            commonLabel_ = "";
            std::cout<<"--> Forced burst mode: Datasets have been ignored. Merge all files."<<std::endl;
        }
    }
    // map assigning a unique key for each dataset  (either jobType or dataset name) and all correponding root files
    std::map<std::string, std::pair<ACDataset*, std::vector<TFile*> > > fileMap;
    std::string allNames = "FORCED";// only filled if forcedMergeAll_
    for (std::vector<std::pair<ACDataset*, TFile*> >::const_iterator pair = allDatasets.begin(); pair != allDatasets.end(); ++pair) {
        //std::cout<<"--> found dataset: "<<dataset->name()<<std::endl;
        std::map<std::string, std::pair<ACDataset*, std::vector<TFile*> > >::iterator entry = fileMap.end();
        if (operationMode_==2) {
            entry = fileMap.find(pair->first->jobType());
        } else if (operationMode_==1) {
            entry = fileMap.find(pair->first->name());
        } else if (operationMode_==3) {
            entry = fileMap.begin();
        }
        if (entry==fileMap.end()) {
            // need a vector here, may not store a TList into a stl container
            std::vector<TFile*> fileList;
            fileList.push_back(pair->second);
            if (operationMode_==2) {
                fileMap.insert(std::make_pair(pair->first->jobType(), std::make_pair(pair->first, fileList)));
            } else if (operationMode_==1) {
                fileMap.insert(std::make_pair(pair->first->name(), std::make_pair(pair->first, fileList)));
            } else if (operationMode_==3) {
                fileMap.insert(std::make_pair("forcedMergeAll", std::make_pair(pair->first, fileList)));
            }
        } else {
            *(entry->second.first) += *(pair->first);//merge datasets, how does this work in forcedMergeAll?
            entry->second.second.push_back(pair->second);//add root file to merge
        }
        if (operationMode_==3) allNames += "_"+pair->first->name();
    }
    if (operationMode_==3) {
        if (fileMap.size()==1) {
            std::map<std::string, std::pair<ACDataset*, std::vector<TFile*> > > newFileMap;
            newFileMap.insert(std::make_pair(allNames, fileMap.begin()->second));
            fileMap = newFileMap;
        }
    }
    mergedDatasets_->clear();
    std::vector<std::string> validStoragePaths;
    for (std::map<std::string, std::pair<ACDataset*, std::vector<TFile*> > >::const_iterator sample = fileMap.begin(); sample != fileMap.end(); ++sample) {
        TFile * file = 0;
        
        std::string currentPath = sample->second.second.front()->GetPath();
        
        if (sample->second.second.size() == 1 && (currentPath.find("merged_"+sample->first+".root")!=std::string::npos || currentPath.find("merged_FORCED_")!=std::string::npos)) {
            // do not merge an input file identical to the merged output
            file = sample->second.second.front();
        } else {
            //merge also single files to have them in the same directory (a copy would be better, I was just too lazy)
            std::string storagePath = outpath_+"/"+commonLabel_;
            bool alreadyTested = false;
            for (std::vector<std::string>::const_iterator validStoragePath = validStoragePaths.begin(); validStoragePath != validStoragePaths.end(); ++validStoragePath) {
                if (*validStoragePath == storagePath) {
                    alreadyTested = true;
                    break;
                }
            }
            if (!alreadyTested) {
                if (!createPlotDir(storagePath)) {
                    std::cout<<"ACHarvester::loadRootFiles:ERROR! Output directory not available!"<<std::endl;
                    throw 404;            
                } else {
                    validStoragePaths.push_back(storagePath);
                }
            }
            std::string fullPath = storagePath+"/merged_"+sample->first+".root";
            file = new TFile(fullPath.c_str(), "RECREATE");
            TList filesOfCommonSample;
            for (std::vector<TFile*>::const_iterator iter = sample->second.second.begin(); iter != sample->second.second.end(); ++iter) {
                filesOfCommonSample.Add(*iter);
            }            
            if (forcedMergeAll_) {
                // choose new name for merged datasets (need to modify the dataset in the ROOT file, not the one in this map)
                std::cout<<"--> Forced to merge datasets. Choose combined name for the merged dataset:"<<std::endl;
                std::cout<<"\tcurrent merged dataset name = "<<sample->second.first->name()<<std::endl;
                std::cout<<"\tType new name (just hit enter to keep the current value): ";
                std::string newValue = readInput();
                if (newValue!="") {
                    forcedMergeName_ = newValue;
                    std::cout<<"\t... replaced by "<<newValue<<std::endl;
                } else std::cout<<"\t... keep it."<<std::endl;
            }
            mergeRootFiles(file, &filesOfCommonSample, 0, 0);
            file->Write();
            // the following two lines are a HOTFIX needed due to a ROOT bug on larger amounts of files
            file->Close();
            if (!forcedMergeAll_) file = TFile::Open(fullPath.c_str(), "READ");//no writing needed any longer, no file needed for forcedMergeAll_, as program will quit
        }
        if (!forcedMergeAll_) {
            if (file) mergedDatasets_->push_back(std::make_pair(sample->second.first, file));
            else {
                std::cout<<"ACHarvester::loadRootFiles:ERROR! Invalid file created!"<<std::endl;
                throw 404;
            }
        }
        if (sample->second.second.size()>1) printf("--> Merged %lu file(s) for %s.\n", sample->second.second.size(), sample->first.c_str());
    }
    if (forcedMergeAll_) {        
        std::cout<<"    Skip further processing in forced merge mode."<<std::endl;
        return false;
    }
    return true;
}
void ACHarvester::burstCompare() {
    if (operationMode_==3) return;
    if (mergedDatasets_->size() != 2) {
        std::cout<<"ACHarvester::burstCompare:INFO! Provide exact two different samples to compare! Skip it."<<std::endl;
        return;
    }
    if (operationMode_==1) printf("--> Start comparison of %s vs. %s.\n", mergedDatasets_->at(0).first->name().c_str(), mergedDatasets_->at(1).first->name().c_str());
    else if (operationMode_==2) printf("--> Start comparison of %s vs. %s.\n", mergedDatasets_->at(0).first->jobType().c_str(), mergedDatasets_->at(1).first->jobType().c_str());

    std::string storagePath = outpath_+"/"+commonLabel_;
    if (!createPlotDir(storagePath)) {
        std::cout<<"ACHarvester::burstCompare:ERROR! Output directory not available!"<<std::endl;
        throw 404;            
    }
    std::string filePath = storagePath+"/comparison_";
    if (operationMode_==1) {
        filePath += mergedDatasets_->at(0).first->name()+"_"+mergedDatasets_->at(1).first->name()+".root";
    } else if (operationMode_==2) {
        filePath += mergedDatasets_->at(0).first->jobType()+"_"+mergedDatasets_->at(1).first->jobType()+".root";
    }
    TFile * target = new TFile(filePath.c_str(), "RECREATE");
    
    TList fileList;
    for (std::vector<std::pair<ACDataset*, TFile*> >::const_iterator sample = mergedDatasets_->begin(); sample != mergedDatasets_->end(); ++sample) {
        fileList.Add(sample->second);
    }

    compareStyle();
    mergeRootFiles(target, &fileList, 0, 1);
    target->Write();
    target->Close();
    printf("    Comparison output at %s.\n", filePath.c_str());
}
void ACHarvester::scanFileForHistograms(const TFile * file) {
    scanRootFile(file, 0);
    //std::cout<<"ACHarvester::scanFileForHistograms done for "<<file->GetName()<<std::endl;
    //file->Close();// User may not want to invalidate the file here
    //delete file;
}

void ACHarvester::mergeHistograms(TH1 * hist1, TH1 * hist2) {
    hist1->Add(hist2);
}
void ACHarvester::compareHistograms(TH1 * hist1, TH1 * hist2, const ACDataset * dataset1, const ACDataset * dataset2, const std::string & storagePath, const std::string & path) {
    if (!createPlotDir(storagePath, path)) {
        std::cout<<"ACHarvester::compareHistograms: ERROR! Failed to create path: "<<storagePath<<"/"<<path<<std::endl;
        throw 404;
    }
    if (hist1->GetEntries() == 0. || hist2->GetEntries() == 0.) {
        std::cout<<"ACHarvester::compareHistograms: WARNING! Empty histogram. Skip "<<path<<"/"<<hist1->GetName()<<std::endl;
        return;
    }
    if (!checkCompatibleBinning(hist1, hist2)) {
        std::cout<<"    Bins mismatch! Skip histogram "<<path<<"/"<<hist1->GetName()<<std::endl;
        return;
    }
//    TH2 * test = dynamic_cast<TH2*>(hist1);
//    if (test) {
//        std::cout<<"    Skip 2d histogram "<<path<<"/"<<hist1->GetName()<<std::endl;
//        return;
//    }
    //std::cout<<"histogram "<<path<<"/"<<hist1->GetName()<<std::endl;
    
    TH1 * h1 = hist1;
    TH1 * h2 = hist2;
    const ACDataset * d1 = dataset1;
    const ACDataset * d2 = dataset2;
    if (!d1->isData() && d2->isData()) {
        // if available: MC solid, data as points
        h1 = hist2;
        h2 = hist1;
        d1 = dataset2;
        d2 = dataset1;
    }
    
    TCanvas	can(h1->GetName(), h1->GetTitle(), 400+4, 500+28);
    TPad mainPad(h1->GetName(),h1->GetName(), 0., 0.25, 1., 1., 0, 0);
    TPad ratioPad("ratio","ratio", 0., 0., 1, .25, 0, 0);
    mainPad.Draw();
    ratioPad.Draw();
    
    int color = 14;

    double integral = h2->Integral(0,h2->GetNbinsX()+1);
    // scale if required
    if (normType_>0) {
        if (!h1->GetSumw2N()) h1->Sumw2();
        if (normType_==1) {
            //scale to common integral
            h1->Scale(1./h1->Integral(0,h1->GetNbinsX()+1));// cannot use GetSumOfWeights as under/overflow would't be considered
            h2->Scale(1./h2->Integral(0,h2->GetNbinsX()+1));
        }
        if (normType_==2) {
            //scale to lumi
            std::cout<<"ACHarvester::compareHistograms: WARNING! Normalization to lumi not implemented yet. Need to reed desired lumi"<<std::endl;
            //h1->Scale(d1->norm(desiredLumi));
            //h2->Scale(d2->norm(desiredLumi));
        }
    }
    
    mainPad.cd();
    h2->UseCurrentStyle();
    h2->SetLineColor(color);//or darker line?
    h2->SetFillColor(color+2);
//    h2->SetMarkerColor(color+2);
//    h2->SetMarkerStyle(kFullSquare);
    h2->SetMarkerSize(0);
    h2->SetFillStyle(1001);//solid
    h1->UseCurrentStyle();
    h1->SetMarkerColor(color);
    if (integral < 100 * h2->GetNbinsX()) {
        h2->DrawCopy("e3");//only draw error contour for reasonably large errors (estimate their average size by unnormalized integral)
        h2->SetFillStyle(0);
        h2->Draw("hist same");
    } else {
        h2->SetFillStyle(0);
        h2->Draw("hist");
    }
    h1->Draw("P0 e1x0 same");
    
    std::vector<TH1*> histVect;
    histVect.push_back(h1);
    histVect.push_back(h2);
    std::vector<std::string> names;
    if (operationMode_==1) {
        names.push_back(d1->name());
        names.push_back(d2->name());
    } else if (operationMode_==2) {
        names.push_back(d1->jobType());
        names.push_back(d2->jobType());    
    }
    addHighLegend(&can, histVect, names);
    addBranding(commonLabel_);
    
    ratioPad.cd();
    TH1 * ratio = (TH1*) h1->Clone("ratio");
    ratio->SetTitle("");
    ratio->GetXaxis()->SetLabelSize(0.);
    ratio->GetXaxis()->SetTitleSize(0.);
    ratio->GetXaxis()->SetTickLength(3*ratio->GetXaxis()->GetTickLength());
    ratio->GetYaxis()->SetTitle("Ratio");
    ratio->GetYaxis()->SetLabelSize(3*ratio->GetYaxis()->GetLabelSize());
    ratio->GetYaxis()->SetTitleOffset(.6);
    ratio->GetYaxis()->SetTitleSize(3*ratio->GetYaxis()->GetTitleSize());
    
    ratio->Divide(h2);
    ratio->SetAxisRange(0., 2., "Y");
    ratio->Draw("P0 e1x0");
    
    int defaultValue = gErrorIgnoreLevel;
	gErrorIgnoreLevel=2000;
    can.Print((storagePath+"/"+path+"/"+can.GetName()+".eps").c_str(),"eps");
    gErrorIgnoreLevel=defaultValue;
    
    delete ratio;//or store ratio into file?
}
bool ACHarvester::modifyDataset(ACDataset * dataset) {
    // do whatever you want to manipulate stored datasets
    std::cout<<"ACHarvester::modifyDataset:"<<std::endl;
    bool madeChanges = false;
    std::cout<<"\tdataset name = "<<dataset->name()<<std::endl;
    std::cout<<"\tType new name (just hit enter to keep the current value): ";
    std::string newValue = readInput();
    if (newValue!="") {
        dataset->setName(newValue);
        std::cout<<"\t... replaced by "<<newValue<<std::endl;
        madeChanges = true;
    } else std::cout<<"\t... keep it."<<std::endl;

    std::cout<<"\tjobType = "<<dataset->jobType()<<std::endl;
    std::cout<<"\tType new jobType (just hit enter to keep the current value): ";
    newValue = readInput();
    if (newValue!="") {
        dataset->setJobType(newValue);
        std::cout<<"\t... replaced by "<<newValue<<std::endl;
        madeChanges = true;
    } else std::cout<<"\t... keep it."<<std::endl;

    return madeChanges;
}
int ACHarvester::mergeRootFiles(TDirectory * target, const TList * sourcelist, int depth, int mode) {
    //std::cout << "Target path: " << target->GetPath() << std::endl;
    if (mode==1) {
        if (sourcelist->GetSize()!=2) {
            std::cout<<"ACHarvester::mergeRootFiles: ERROR! To compare histograms specify exact 2 files, not "<<sourcelist->GetSize()<<std::endl;
            throw 404;
        }
    }
    depth++;
//    printf("--> Merge layer %i\r", depth);
//    fflush(stdout);

    
    Bool_t status = kTRUE;
    
    // Get the dir name
    std::string fullPath = target->GetPath();
    std::string::size_type pos = fullPath.find(":");
    TString path(fullPath.substr(pos+2));
    std::string storagePath = fullPath.substr(0, pos-5);//strip .root
    //std::cout << "storagePath, path: " << storagePath <<","<< path << std::endl;

    
    Int_t nguess = sourcelist->GetSize()+1000;
    THashList allNames(nguess);
    allNames.SetOwner(kTRUE);//If hashlist is the owner of its contents, these objects will be deleted whenever the collection itself is deleted

    // performence increase
    ((THashList*)target->GetList())->Rehash(nguess);
    ((THashList*)target->GetListOfKeys())->Rehash(nguess);
    
    
    TFile      * current_file = (TFile*)sourcelist->First();
    TDirectory * current_sourcedir = current_file->GetDirectory(path);
    
    while (current_file || current_sourcedir) {
        // When current_sourcedir != 0 and current_file == 0 we are going over the target
        // for an incremental merge.
        if (current_sourcedir && (current_file == 0 || current_sourcedir != target)) {
            
            // loop over all keys in this directory
            TIter nextkey( current_sourcedir->GetListOfKeys() );
            TKey *key;
            TString oldkeyname;
            
            while ( (key = (TKey*)nextkey())) {
                
                // Keep only the highest cycle number for each key.  They are stored in the (hash) list
                // consecutively and in decreasing order of cycles, so we can continue until the name
                // changes.
                if (oldkeyname == key->GetName()) continue;
                // Read in but do not copy directly the processIds.
                if (strcmp(key->GetClassName(),"TProcessID") == 0) { key->ReadObj(); continue;}
                // If we have already seen this object [name], we already processed
                // the whole list of files for this objects and we can just skip it
                // and any related cycles.
                if (allNames.FindObject(key->GetName())) {
                    oldkeyname = key->GetName();
                    continue;
                }
                
                TClass *cl = TClass::GetClass(key->GetClassName());
                if (!cl || !cl->InheritsFrom(TObject::Class())) {
                    Info("MergeRecursive", "cannot merge object type, name: %s title: %s",
                         key->GetName(), key->GetTitle());
                    continue;
                }
                allNames.Add(new TObjString(key->GetName()));
                
                // read object from first source file
                TObject *obj;
                obj = key->ReadObj();
                if (!obj) {
                    Info("MergeRecursive", "could not read object for key {%s, %s}",
                         key->GetName(), key->GetTitle());
                    continue;
                }
                
                Bool_t canBeMerged = kTRUE;
                
                if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
                    // it's a subdirectory
                    //std::cout << "Found subdirectory " << obj->GetName() << "at depth "<< depth << std::endl;
                    target->cd();
                    TDirectory *newdir;
                    newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
                    
                    // newdir is now the starting point of another round of merging
                    // newdir still knows its depth within the target file via
                    // GetPath(), so we can still figure out where we are in the recursion
                    //status = MergeRecursive(newdir, sourcelist, type);
                    status = mergeRootFiles(newdir, sourcelist, depth, mode);
                    if (!status) return status;
                } else if (obj->IsA()->InheritsFrom("TH1") || obj->IsA()->InheritsFrom("ACDataset") ) {
                    TList inputs;//store link to input objects for easy deletion

                    ACDataset * dataset = 0;
                    TH1 * h1 = 0;
                    
                    if (obj->IsA()->InheritsFrom("TH1") ) {
                        h1 = (TH1*)obj;
                        //std::cout << "Work on histogram " << obj->GetName() << std::endl;
                    } else if (obj->IsA()->InheritsFrom("ACDataset") ) {
                        dataset = (ACDataset*)obj;
                        std::cout<<"--> found first dataset: "<<dataset->name()<<std::endl;
                    }

                    // Loop over all source files and merge same-name object
                    TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
                    if (nextsource == 0) {
                        // There is only one file in the list
                    } else {
                        unsigned int sampleID = 0;
                        bool inconsistency = false;
                        bool binInconsistency = false;
                        do {
                            sampleID++;
                            // make sure we are at the correct directory level by cd'ing to path
                            TDirectory *ndir = nextsource->GetDirectory(path);
                            if (ndir) {
                                ndir->cd();
                                TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
                                if (key2) {
                                    TObject *hobj = key2->ReadObj();
                                    if (!hobj) {
                                        Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
                                             key->GetName(), key->GetTitle(), nextsource->GetName());
                                        nextsource = (TFile*)sourcelist->After(nextsource);
                                        continue;
                                    }
                                    // Set ownership for collections
                                    if (hobj->InheritsFrom(TCollection::Class())) {
                                        ((TCollection*)hobj)->SetOwner();
                                    }
                                    hobj->ResetBit(kMustCleanup);
                                    inputs.Add(hobj);
                                    
                                    // now do indiviual merging operations
                                    if (obj->IsA()->InheritsFrom("TH1") ) {
                                        TH1 *h2 = (TH2*)hobj;
                                        if (mode==0) {
                                            if (checkCompatibleBinning(h1, h2)) mergeHistograms(h1, h2);
                                            else {
                                                if (!recoverBinInconsistency(h1, h2)) {
                                                    if (!binInconsistency) {
                                                        // report inconsistencies once
                                                        std::cout<<"    Bins mismatch! Skip histogram "<<path<<"/"<<h2->GetName()<<". No further warnings."<<std::endl;
                                                        binInconsistency = true;
                                                    }
                                                }
                                            }
                                        } else compareHistograms(h1, h2, mergedDatasets_->at(0).first, mergedDatasets_->at(sampleID).first, storagePath, path.Data());
                                    }
                                    if (obj->IsA()->InheritsFrom("ACDataset") ) {
                                        if (mode==0) {
                                            ACDataset * dataset2 = (ACDataset*)hobj;
                                            //std::cout<<"--> found further dataset: "<<dataset2->name()<<std::endl;
                                            if (!forcedMergeAll_) {
                                                if (dataset->name()!=dataset2->name()) {
                                                    std::cout<<"ACHarvester::mergeRootFiles:ERROR! Mismatch in ACDatasets names: "<<dataset->name()<<"!="<<dataset2->name()<<std::endl;
                                                    throw 404;
                                                }
                                                if (dataset->jobType()!=dataset2->jobType()) {
                                                    std::cout<<"ACHarvester::mergeRootFiles:ERROR! Mismatch in ACDataset jobTypes: "<<dataset->jobType()<<"!="<<dataset2->jobType()<<std::endl;
                                                    throw 404;
                                                }
                                            }
                                            combineDatasets(dataset, dataset2);
                                        }
                                    }
                                    // Remove all objects from the list AND delete all heap based objects
                                    inputs.Delete();
                                } else {
                                    if (!inconsistency) {
                                        // report inconsistencies once
                                        std::cout<<"    Skip missing object "<<obj->GetName()<<". No further warnings."<<std::endl;
                                        inconsistency = true;
                                    }
                                }
                            } //std::cout<<"    ignore missing dir "<<path<<" for obj "<<obj->GetName()<<" in file "<<nextsource->GetPath()<<std::endl;
                            nextsource = (TFile*)sourcelist->After( nextsource );
                        } while (nextsource);
                    }
                } else {
                    // Object is of no type that we can merge
                    std::cout << "ACHarvester::mergeRootFiles: WARNING! Cannot merge unknown object type, name: "<<obj->GetName() << " title: " << obj->GetTitle() << std::endl;

                    Bool_t warned = kFALSE;
                    canBeMerged = kFALSE;
                    
                    // Loop over all source files and write similar objects directly to the output file
                    TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
                    while (nextsource) {
                        // make sure we are at the correct directory level by cd'ing to path
                        TDirectory *ndir = nextsource->GetDirectory(path);
                        if (ndir) {
                            ndir->cd();
                            TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
                            if (key2) {
                                if (!warned) {
                                    Warning("MergeRecursive", "cannot merge object type (n:'%s', t:'%s') - "
                                            "Merge(TCollection *) not implemented",
                                            obj->GetName(), obj->GetTitle());
                                    warned = kTRUE;
                                }
                                TObject *nobj = key2->ReadObj();
                                if (!nobj) {
                                    Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
                                         key->GetName(), key->GetTitle(), nextsource->GetName());
                                    nextsource = (TFile*)sourcelist->After(nextsource);
                                    continue;
                                }
                                nobj->ResetBit(kMustCleanup);
                                if (target->WriteTObject(nobj, key2->GetName(), "SingleKey") <= 0) {
                                    Warning("MergeRecursive", "problems copying object (n:'%s', t:'%s') to output file ",
                                            obj->GetName(), obj->GetTitle());
                                    status = kFALSE;
                                }
                                delete nobj;
                            }
                        }
                        nextsource = (TFile*)sourcelist->After( nextsource );
                    }
                }
                
                // now write the merged histogram (which is "in" obj) to the target file
                // note that this will just store obj in the current directory level,
                // which is not persistent until the complete directory itself is stored
                // by "target->SaveSelf()" below
                target->cd();
                
                oldkeyname = key->GetName();
                //!!if the object is a tree, it is stored in globChain...
                if(obj->IsA()->InheritsFrom( TDirectory::Class() )) {
                    // Do not delete the directory if it is part of the output
                    // and we are in incremental mode (because it will be reuse
                    // and has not been written to disk (for performance reason).
                    // coverity[var_deref_model] the IsA()->InheritsFrom guarantees that the dynamic_cast will succeed. 
                    if (dynamic_cast<TDirectory*>(obj)->GetFile() != target) {
                        delete obj;
                    }
                } else if (obj->IsA()->InheritsFrom( TCollection::Class() )) {
                    // Don't overwrite, if the object were not merged.
                    if ( obj->Write( oldkeyname, canBeMerged ? TObject::kSingleKey | TObject::kOverwrite : TObject::kSingleKey) <= 0 ) {
                        status = kFALSE;
                    }
                    ((TCollection*)obj)->SetOwner();
                    delete obj;
                } else {
                    // Don't overwrite, if the object were not merged.
                    if ( obj->Write( oldkeyname, canBeMerged ? TObject::kOverwrite : 0) <= 0) {
                        status = kFALSE;
                    }
                    delete obj;
                }
            } // while ( ( TKey *key = (TKey*)nextkey() ) )
        }
        current_file = current_file ? (TFile*)sourcelist->After(current_file) : (TFile*)sourcelist->First();
        if (current_file) {
            current_sourcedir = current_file->GetDirectory(path);
        } else {
            current_sourcedir = 0;
        }
    }
    // save modifications to the target directory.

    // In case of incremental build, we will call Write on the top directory/file, so we do not need
    // to call SaveSelf explicilty.
    target->SaveSelf(kTRUE);
    
    return status;
}
void ACHarvester::scanRootFile(const TDirectory * target, int depth) {
    //std::cout << "Target path: " << target->GetPath() << std::endl;
    depth++;
    //    printf("--> Merge layer %i\r", depth);
    //    fflush(stdout);
    
    //gain time, do not add the objects in the list in memory
    Bool_t status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);//It is now your responsibility to delete the following histograms once you have finished with them.
    
    // loop over all keys in this directory
    TIter nextkey(target->GetListOfKeys());
    TKey *key, *oldkey=0;
    while ((key = (TKey*)nextkey())) {
        //std::cout << "move to next oject at " << key << std::endl;
        //keep only the highest cycle number for each key
        if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
        // read object from first source file
        TObject *obj = key->ReadObj();
        //std::cout << "Reading oject " << obj->GetName() << std::endl;
        
        if (obj->IsA()->InheritsFrom("TDirectory")) {
            //std::cout << "Found subdirectory " << obj->GetName() << " at depth "<< depth << std::endl;
            //TDirectory * newdir = target->GetDirectory(obj->GetName());
            scanRootFile((TDirectory*)obj, depth);
        } else if (obj->IsA()->InheritsFrom("TH1")) {
            TH1 *h1 = (TH1*)obj;
            //std::cout << "Process histogram " << h1->GetName() << " with nentries " << h1->GetEntries() << std::endl;
            std::string fullPath = target->GetPath();
            std::string::size_type pos = fullPath.find(":");
            TString path(fullPath.substr(pos+2));
            processHistogram(h1, path.Data());
        }
        oldkey = key;        
    }
    TH1::AddDirectory(status);
}
void ACHarvester::combineDatasets(ACDataset * d1, const ACDataset * d2) {
    //std::cout<<"merge dataset: "<<*d1<<" and "<<*d2<<std::endl;
    //std::cout<<"merge dataset: "<<d1->name()<<" and "<<d2->name()<<std::endl;
    *d1 += *d2;
    if (forcedMergeAll_ && forcedMergeName_ != "") {
        d1->setName(forcedMergeName_);
    }
}
bool ACHarvester::recoverBinInconsistency(TH1 * h1, TH1 * h2) const {
    bool solved = false;
    //user-defined solution expected here
    
    return solved;
}

// helper functions needed in compareHistograms
void ACHarvester::scanDirectory(const std::string & directory, std::vector<std::string> & fileNames, int level) {
    if (level==0) {
        printf("--> Scan directory %s ...\n", directory.c_str());
    }
    if (boost::filesystem::is_directory(directory)) {
        level++;
        for (boost::filesystem::directory_iterator itr(directory); itr!=boost::filesystem::directory_iterator(); ++itr) {
            //std::cout<<itr->path().string()<<std::endl;// display full path
            if (boost::filesystem::is_directory(itr->path().string())) {
                printf("%*s cd %s\n", 4*level-1, "", itr->path().filename().c_str());// display filename only
                scanDirectory(itr->path().string(), fileNames, level);
            } else if (is_regular_file(itr->status())) {
                if (itr->path().extension()==".root") {
                    printf("%*s adding %s ", 4*level-1, "", itr->path().filename().c_str());// display filename only
                    std::cout<<"["<<file_size(itr->path())/1000<<" kB]"<<std::endl;
                    fileNames.push_back(itr->path().string());
                }
            }
        }
    } else std::cout << (boost::filesystem::exists(directory) ? "    Found: " : "    Not found: ") << directory << '\n';
}
bool ACHarvester::createPlotDir(const std::string & storagePath, const std::string & path) const {
    struct stat st;
    if (stat(storagePath.c_str(), &st)!=0) {
        std::cout<<"ACHarvester::createPlotDir: the following path does not exist\n\t"<<storagePath<<std::endl;
        std::cout<<"\tCreate it? (y/n)\t";
        std::string create = readInput();
        if (create=="y") {
            if (mkdir(storagePath.c_str(), 0777) == 0) std::cout<<"\t"<<storagePath<<" was created."<<std::endl;
            else{
                std::cout<<"\tFailed to create "<<storagePath<<"."<<std::endl;
                return false;
            }
        } else return false;
    } else {
        std::cout<<"ACHarvester::createPlotDir: the following path already exists\n\t"<<storagePath<<std::endl;
        std::cout<<"\tKeep it and write into the existing folder? (y/n)\t";
        std::string create = readInput();
        if (create!="y") return false;
    }
    
    if (path=="") return true;

	//create subdirs (including /'s in canName)
    std::vector<std::string> folders;
    tokenizePath(path, folders, "/");
	std::string tmp = "";
    for (std::vector<std::string>::const_iterator folder = folders.begin(); folder != folders.end(); ++folder) {
        tmp += "/" + *folder;
		if (stat((storagePath+"/"+tmp).c_str(), &st)==0) {//already exists, so continue
			continue;
		}
		//std::cout<<"Plot::setupDirectories:Info: Path "<<path_+tmp<<" does not exist ";
		if (mkdir((storagePath+"/"+tmp).c_str(), 0777) != 0) {
			std::cout<<"ACHarvester::createPlotDir: Path "<<storagePath+tmp<<" does not exist and failed to be created.\n";
			break;//break if failed to create
		}// else std::cout<<"but was created.\n";
    }
    
    return true;
}
bool ACHarvester::prepareModificationOfDataset(TFile * file, TKey * key, ACDataset * dataset) {
    bool modified = false;
    if (!modifyInputFiles_) return modified;

    if (modifyDataset(dataset)) {
        dataset->Write("", TObject::kOverwrite);
        file->Write();
        modified = true;
    }
    return modified;
}
bool ACHarvester::checkCompatibleBinning(const TH1 * h1, const TH1 * h2) const {
    if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY() || h1->GetNbinsZ() != h2->GetNbinsZ()) {
        return false;
    }
    if (!checkAxisLimits(h1->GetXaxis(), h2->GetXaxis()) || !checkAxisLimits(h1->GetYaxis(), h2->GetYaxis()) || !checkAxisLimits(h1->GetZaxis(), h2->GetZaxis())) {
        return false;
    }
    
    return true;
}
bool ACHarvester::checkAxisLimits(const TAxis *a1, const TAxis *a2) const {
    if (! TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) || ! TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12)) {
        return false;
    }
    return true;
}
void ACHarvester::addHighLegend(TVirtualPad *can, const std::vector<TH1*> & histVect, const std::vector<std::string> & names) {
	unsigned int columns = histVect.size();
	if(columns==0) return;
    if (columns!=names.size()) {
        std::cout<<"ACHarvester::addHighLegend:ERROR! Size mismatch between histograms and labels. Skip legend."<<std::endl;
        return;
    }
	float wholeWidth = 1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin();
	float height = 0.05;
	float x1 = gStyle->GetPadLeftMargin()-.0075;//offset due to padding within TLegend
	float y1 = 1-gStyle->GetPadTopMargin();
	float x2 = wholeWidth+gStyle->GetPadLeftMargin();
	float y2 = 1-gStyle->GetPadTopMargin()+height;
	TLegend *legend = new TLegend(x1, y1, x2, y2,"","BRNDC");
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	//		legend->SetTextFont(62);
	legend->SetNColumns(columns);
	legend->SetMargin(height/wholeWidth*columns);//create quadratic symbols
	legend->SetEntrySeparation(0.00);
	legend->SetColumnSeparation(0.00);
	legend->SetTextAlign(12);//3*10=right,3*1=top
	unsigned int index = 0;
	TLegendEntry * entry;
	for(std::vector<TH1*>::const_iterator iter=histVect.begin(); iter!=histVect.end(); ++iter, index++) {
		std::string label = names.at(index);
        std::string drawOption = "F";
        if(std::string((*iter)->GetDrawOption()).find("P",0,1)!=std::string::npos) drawOption = "P";
        entry = legend->AddEntry(*iter, label.c_str(), drawOption.c_str());
		legend->SetTextSize(gStyle->GetTitleSize());//same size for labels and legend
	}
	legend->Draw();
	//		delete legend;//has to survive
}
void ACHarvester::addBranding(const std::string & branding) {
    if (branding=="") return;
    
	TLatex *text = new TLatex;
	text->SetNDC();
	text->SetTextAlign(13);
	text->SetTextSize(.025);
	text->SetTextColor(13);
	
    text->DrawLatex(.005, .995, branding.c_str());
	
	delete text;
}
void ACHarvester::compareStyle(){
	int fontID = 42;//4=helvetica-medium-r-normal + precision 2, ROOT default is 62
    
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	gStyle->SetLineScalePS(2);
    gStyle->SetPaperSize(20,25);
	
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.07);
	gStyle->SetPadBorderSize(0);
	
	gStyle->SetTickLength(0.03, "xyz");
	gStyle->SetLabelSize(0.035,"xyz");
	gStyle->SetTitleSize(0.0379,"xyz");
	gStyle->SetTitleOffset(1.85,"yz");
	gStyle->SetTitleOffset(1.2,"x");
	gStyle->SetLabelOffset(0.008,"yz");
	gStyle->SetLabelOffset(0.01,"x");
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLabelFont(fontID, "xyz");
	
	gStyle->SetStatX(1-gStyle->GetPadRightMargin());
	gStyle->SetStatY(1-gStyle->GetPadTopMargin());
	gStyle->SetStatW(0.2);
	gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(fontID);
	
	gStyle->SetTitleAlign(13);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleX(gStyle->GetPadLeftMargin());
	gStyle->SetTitleY(.995);
	gStyle->SetTitleW(1-gStyle->GetPadRightMargin()-gStyle->GetPadLeftMargin());
	gStyle->SetTitleH(.04);
	gStyle->SetTitleFont(fontID, "title");
	gStyle->SetTitleFont(fontID, "xyz");
	
	gStyle->SetTitleFontSize(0.01);
	gStyle->SetTextAlign(12);
	gStyle->SetTextFont(fontID);
	
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	
	gStyle->SetCanvasColor(0);
	gStyle->SetPadColor(0);
	
	gStyle->SetHistFillStyle(0);
	gStyle->SetHistLineColor(1);
	gStyle->SetHistLineWidth(1);
	
	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
	gStyle->SetGridColor(17);
	gStyle->SetGridStyle(2);
	gStyle->SetGridWidth(1);
	
	gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(.5);
}
void ACHarvester::tokenizePath(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) const {
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
const std::string ACHarvester::readInput() const {
    std::string result = "";
    // clear all recent status flags to enable further reading
    std::cin.clear();
    // Ignore to the end of file
    //std::cin.ignore(INT_MAX); // or std::cin.ignore(std::numeric_limits<std::streamsize>::max());
    // Ignore to the end of line
    //std::cin.ignore(INT_MAX, '\n'); //or std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // do not skip leading whitespace
    std::cin.unsetf(std::ios_base::skipws);
 
    if (!(std::cin >> result)) {
		std::cin.clear();// clear again in case of errors
	}
	while (std::cin.get() != '\n');// read until newline
    //std::cout<<"ACHarvester::readInput: result is "<<result<<std::endl;
    
    return result;
}
