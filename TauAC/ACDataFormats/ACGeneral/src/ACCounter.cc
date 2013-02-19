#include "../interface/ACCounter.h"

ClassImp(ACCounter)

ACCounter::ACCounter() {
    
}
ACCounter::~ACCounter() {
    
}
void ACCounter::add(const std::string & name, double start) {
    if (!counts_.insert(make_pair(name, start)).second) {
        std::cerr << "Error in ACCounter::add: ACCounter '" << name << "' already exists!" << std::endl;
    }
}
void ACCounter::set(const std::string & name, double value) {
    std::map<std::string,double>::iterator iter = counts_.find(name);
    if (iter != counts_.end()) {
        iter->second = value;
    } else {
        std::cerr << "Error in ACCounter::get: ACCounter '" << name << "' does not exists!" << std::endl;
    }
}
double ACCounter::get(const std::string & name) const {
    std::map<std::string,double>::const_iterator iter = counts_.find(name);
    if (iter != counts_.end()) {
        return iter->second;
    } else {
        std::cerr << "Error in ACCounter::get: ACCounter '" << name << "' does not exists!" << std::endl;
        return -1.0;
    }
    return -1.0;
}
int ACCounter::getInt(const std::string & name) const { return int(get(name)); }
void ACCounter::count(const std::string & name, double amount) {
    std::map<std::string,double>::iterator iter = counts_.find(name);
    if (iter != counts_.end()) {
        (iter->second) += amount;
    } else {
        std::cerr << "Error in ACCounter::count: ACCounter '" << name << "' does not exists!" << std::endl;
    }
}
double ACCounter::ratio(const std::string & numerator, const std::string & denominator, int precision) const {
    double eff = -1.0;
    std::map<std::string,double>::const_iterator iternum = counts_.find(numerator);
    std::map<std::string,double>::const_iterator iterdenom = counts_.find(denominator);
    if (iternum != counts_.end() && iterdenom != counts_.end()) {
        if ((iterdenom->second) > 0 && (iternum->second) >= 0) {
            eff = iternum->second/iterdenom->second;
        }
    } else {
        std::cerr << "Error in ACCounter::ratio: ACCounter '" << numerator << "' or '" << denominator << "' does not exists!" << std::endl;
    }
    
    if (precision == -1000) return eff;
    return round(eff, precision);
}
ACCounter & ACCounter::operator+=(const ACCounter & rhs) {
    for (std::map<std::string, double>::const_iterator iter = rhs.counts_.begin(); iter != rhs.counts_.end(); ++iter) {
        std::map<std::string,double>::iterator found = counts_.find(iter->first);
        if (found != counts_.end()) {
            found->second += iter->second;
        } else {
            counts_.insert(*iter);
        }
    }
    return *this;
}
std::ostream & operator<<(std::ostream & stream, const ACCounter & obj) {
    for (std::map<std::string, double>::const_iterator iter = obj.counts_.begin(); iter != obj.counts_.end(); ++iter) {
        stream << iter->first << " --> " << iter->second << std::endl;
    }        
    return stream;
}
