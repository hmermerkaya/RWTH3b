#ifndef ACCommon_h
#define ACCommon_h
/**
 Base class containing common functions.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include <memory>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <cmath>

#include "TObject.h"

class ACCommon : public TObject {
public:
    /// virtual classname
    virtual std::string classname() const = 0;

    ClassDef(ACCommon,1);

protected:
    /// convert any value to a string
    template <class T> std::string any2string(T i);
    /// split a string into tokens separated by a specific delimiter
    void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
    /// round value to a given precission
    double round(double value, int power = 0);
    /// delete pointers in a vector
    template <class T> void deleteVectorOfPointers(T * inVectorOfPointers);
};

template <class T> std::string ACCommon::any2string(T i) {
    std::ostringstream buffer;
    buffer << i;
    return buffer.str();
}
template <class T> void ACCommon::deleteVectorOfPointers(T * inVectorOfPointers) {
    if (!inVectorOfPointers) {
        return;
    }
    typename T::iterator i;
    for (i = inVectorOfPointers->begin(); i < inVectorOfPointers->end(); ++i) {
        if (!*i) continue;
        delete * i;
        *i = 0;
    }
}

#endif
