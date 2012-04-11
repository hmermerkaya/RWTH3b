#ifndef ACDataset_h
#define ACDataset_h
/**
 Container class to store all information about one dataset.
 
 @author Lars Perchalla & Philip Sauerland
 @date 2012
 */

#include <string>
#include "ACCommon.h"
#include "ACCounter.h"

class ACDataset : public ACCommon {
public:
    virtual std::string classname() const { return "ACDataset"; }
	ACDataset();
	ACDataset(const std::string & name, const std::string & alias, long long int initialEvents, std::vector<std::string> storagePaths, const std::string & jobType, double xsection, float luminosity=-1.0, bool isBackground = false, bool isData = false, double scale = 1.0);
	~ACDataset();
    
	const std::string & name() const { return name_; }
    void setName(const std::string & name) { name_ = name;};
	const std::string & alias() const { return alias_; }
    void setAlias(const std::string & alias) { alias_ = alias;}
	long long int initialEvents() const { return initialEvents_; }
    void setInitialEvents(long long int initialEvents) { initialEvents_ = initialEvents; }
	double xsection() const;
    void setXsection(double xsection) { xsection_ = xsection; }
	float luminosity() const;
    void setLuminosity(float luminosity) { luminosity_ = luminosity; }
	const std::vector<std::string> & paths() const { return storagePaths_; }
    void setStoragePaths(const std::vector<std::string> & storagePaths) { storagePaths_ = storagePaths; }
    const std::vector<int> & runIDs() const { return runIDs_; }
    void setRunIDs(std::vector<int> runIDs) { runIDs_ = runIDs; }
	const std::vector<std::string> & triggerPaths() const { return triggerPaths_; }
    void setTriggerPaths(const std::vector<std::string> & triggerPaths) { triggerPaths_ = triggerPaths; }
	const std::string & jobType() const { return jobType_; }
    void setJobType(const std::string & jobType) { jobType_ = jobType;};
	double scale() const { return scale_; }
    void setScale(const double & scale) { scale_ = scale; }
	bool isBackground() const { return isBackground_; }
    void setIsBackground(bool isBackground) { isBackground_ = isBackground; }
	bool isData() const { return isData_; }
    void setIsData(bool isData) { isData_ = isData; }
    ACCounter * counters() { return counters_; }

    /// obtain normalization factor
	double norm(double desiredLumi) const;
    /// Combine two ACDatasets
    ACDataset & operator+=(const ACDataset & rhs);
    /// overloaded stream operator
    friend std::ostream & operator<<(std::ostream & stream, const ACDataset & obj);

    ClassDef(ACDataset,1);

protected:
	/// unique key of the dataset
	std::string name_;
    /// alias name
    std::string alias_;
	/// number of generated events per sample (ignoring weights)
	long long int initialEvents_;
	/// paths of corresponding filenames for the given dataset (more than one possible due to splitting)
	std::vector<std::string> storagePaths_;
    /// list of processed run ID's
	std::vector<int> runIDs_;
    /// list of processed trigger paths (in case a trigger pattern is selected)
	std::vector<std::string> triggerPaths_;
    /// type of the given job
	std::string jobType_;
	/// cross section in pb, see eventRates.nb or dataSets.numbers. only filled in case of MC.
	double xsection_;
	/// the real indiviual luminosity in pb. only filled in case of data.
	float luminosity_;
	/// true if set is regarded as background
	bool isBackground_;
	/// true if set is real data (and not MC)
	bool isData_;
	/// arbitrary event scale (e.g. maxEvents/nentries, systematics, ...)
	double scale_;
    /// counters
    ACCounter * counters_;
};

typedef ACDataset* pACDataset;

#endif
