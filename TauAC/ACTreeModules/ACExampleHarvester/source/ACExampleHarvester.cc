#include "../interface/ACExampleHarvester.h"

ACExampleHarvester::ACExampleHarvester(const std::string & outpath, bool modifyInputFiles):
ACHarvester(outpath, modifyInputFiles) {
}
ACExampleHarvester::~ACExampleHarvester() {
}

int main(int argc, char* argv[]) {
    TStopwatch timer;
    timer.Start();
    
    if (argc < 2) {
        printf("Usage: %s <whitespace separated list of files/directory containing files>\n", argv[0]);
        return 1;
    }
    //may contain different types of datasets
    std::vector<std::string> allFilenames;
    std::vector<std::string> directories;
    for (int i=1; i!=argc; i++) {
        std::string::size_type pos = std::string(argv[i]).find(".root");
        if (pos!=std::string::npos) {
            if (directories.size()!=0) {
                printf("Usage: %s <whitespace separated list of files/directory containing files>\n both a path and root files are not accepted.\n", argv[0]);
                return 1;
            }
            allFilenames.push_back(argv[i]);
        } else {
            // expect directory now
            if (allFilenames.size()!=0) {
                printf("Usage: %s <whitespace separated list of files/directory containing files>\n both a path and root files are not accepted.\n", argv[0]);
                return 1;
            }
            directories.push_back(argv[i]);
        }
    }
    
    ACExampleHarvester harvester( "/user/perchalla/plots/analysis/CMSSW_4_4_2/Harvester/", false);
    if (allFilenames.size()==0) harvester.scanDirectories(directories, allFilenames);
    harvester.loadRootFiles(allFilenames);
    // work on different datasets
    harvester.burstCompare();
    
    timer.Stop();
    std::cout<<"ACHarvester finished (in "<<timer.CpuTime()<<"s)."<<std::endl;
    
    return 0;
}
