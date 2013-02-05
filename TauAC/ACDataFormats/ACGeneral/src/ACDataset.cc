#include "../interface/ACDataset.h"

ClassImp(ACDataset)

ACDataset::ACDataset() {
	name_ = "";
	alias_ = "";
	initialEvents_ = 0;
	storagePaths_.clear();
    runIDs_.clear();
    triggerPaths_.clear();
	scale_ = 1.0;
	xsection_ = -1.;
	luminosity_ = -1.;
	isBackground_ = false;
	isData_ = false;
    counters_ = new ACCounter();
}
ACDataset::ACDataset(const std::string & name, const std::string & alias, long long int initialEvents, std::vector<std::string> storagePaths, const std::string & jobType, double xsection, float luminosity, bool isBackground, bool isData, double scale):
name_(name), alias_(alias), initialEvents_(initialEvents), storagePaths_(storagePaths), jobType_(jobType), xsection_(xsection), luminosity_(luminosity), isBackground_(isBackground), isData_(isData), scale_(scale) {
    counters_ = new ACCounter();
    runIDs_.clear();
    triggerPaths_.clear();
}
ACDataset::~ACDataset() {
    //delete counters_;
}

double ACDataset::xsection() const {
    if (isData_ && xsection_ == -1.) printf("ACDataset::luminosity: Warning! Your are accessing a value not stored for data.\n");
    return xsection_;
}
float ACDataset::luminosity() const {
    if (!isData_ && luminosity_==-1.) printf("ACDataset::luminosity: Warning! Your are accessing a value not stored for MC.\n");
    return luminosity_;
}
double ACDataset::norm(double desiredLumi) const {
	double ratio = 0.0;
    if (isData_) {
        if (luminosity() > 0.0) ratio = scale()*desiredLumi/luminosity();
    } else {
        if (initialEvents() > 0) ratio = scale()*xsection()*desiredLumi/initialEvents();
        else if (luminosity() > 0.0) ratio = scale()*desiredLumi/luminosity(); // some MC might carry lumi info instead of events (e.g. QCD derived from data)
    }
	
	return ratio;
}
ACDataset & ACDataset::operator+=(const ACDataset & rhs) {
    storagePaths_.insert(storagePaths_.end(), rhs.storagePaths_.begin(), rhs.storagePaths_.end());
    runIDs_.insert(runIDs_.end(), rhs.runIDs_.begin(), rhs.runIDs_.end());
    (*counters_)+=(*rhs.counters_);
    // skip duplicates (should move to map here)
    for (std::vector<std::string>::const_iterator newpath = rhs.triggerPaths_.begin(); newpath != rhs.triggerPaths_.end(); ++newpath) {
        bool duplicate = false;
        for (std::vector<std::string>::const_iterator path = triggerPaths_.begin(); path != triggerPaths_.end(); ++path) {
            if (*path == *newpath) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) triggerPaths_.insert(triggerPaths_.end(), *newpath);
    }
    
    return *this;
}
std::ostream & operator<<(std::ostream & stream, const ACDataset & obj) {
    stream << "dataset name: " << obj.name_ << std::endl;
    stream << "alias: " << obj.alias_ << std::endl;
    stream << "generated events: " << obj.initialEvents_ << std::endl;
    stream << "job type: " << obj.jobType_ << std::endl;
    stream << "cross section: " << obj.xsection_ << " pb" << std::endl;
    stream << "luminosity: " << obj.luminosity_ << " pb^-1" << std::endl;
    stream << "is background: " << obj.isBackground_ << std::endl;
    stream << "is data: " << obj.isData_ << std::endl;
    stream << "scale: " << obj.scale_ << std::endl;
    stream << "counters: " << std::endl;
    stream << *obj.counters_ << std::endl;
    stream << "storage paths: " << std::endl;
    for (std::vector<std::string>::const_iterator iter = obj.storagePaths_.begin(); iter != obj.storagePaths_.end(); ++iter) {
        stream << " " << *iter << std::endl;
    }
    return stream;
}
