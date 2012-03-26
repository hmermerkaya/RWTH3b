#ifndef ACCounter_h
#define ACCounter_h
/**
 Class to count an arbitrary number of things.
 
 @author Lars Perchalla & Philip Sauerland
 @date 2012
 */

#include <map>
#include "ACCommon.h"

class ACCounter : public ACCommon {
public:
    ACCounter();
    virtual ~ACCounter();
    virtual std::string classname() const { return "ACCounter"; }
    void add(const std::string & name, double start=0.0);
    void set(const std::string & name, double value);
    double get(const std::string & name) const;
    int getInt(const std::string & name) const;
    void count(const std::string & name, double amount=1.0);
    double ratio(const std::string & numerator, const std::string & denominator) const;
    const std::map<std::string, double> & counts() const {return counts_;}
    /// Combine two ACCounters
    ACCounter & operator+=(const ACCounter & rhs);
    /// overloaded stream operator
    friend std::ostream & operator<<(std::ostream & stream, const ACCounter & obj);
    
    ClassDef(ACCounter,1);

private:
    std::map<std::string, double> counts_;
};

typedef ACCounter* pACCounter;

#endif
