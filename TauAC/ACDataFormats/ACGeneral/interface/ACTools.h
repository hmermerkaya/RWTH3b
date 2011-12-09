#ifndef ACTools_h
#define ACTools_h
/**
 This class provides the helper functions of ACCommon, whenever no instance of a derived object is available. It does not store any members and shall not have additional funtionality. Only static templates, which need to be call as object function instead of a member function, are collected here.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCommon.h"

class ACTools : public ACCommon {
public:
    virtual std::string classname() const { return "ACTools"; }
    ACTools() {};
    virtual ~ACTools() {};

    /// compare pt. function to be used in std::sort.
    template <class T> static bool cmpPt(const T &a, const T &b) {
        return a->pt() > b->pt();
    }
    /// compare pz. function to be used in std::sort.
    template <class T> static bool cmpPz(const T &a, const T &b) {
        return a->pz() > b->pz();
    }

    ClassDef(ACTools,1);
};

#endif
