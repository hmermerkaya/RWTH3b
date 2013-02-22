#ifndef ACProvenance_h
#define ACProvenance_h
/**
 This class contains the branches of the main tree.
 
 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../../ACDataFormats/ACGeneral/interface/ACCounter.h"

class ACTreeReader;

class ACProvenance {
public:
    ACProvenance();

    /// return content merged over all chain elements
    const ACCounter & counter() const { return mergedCounter_; }
    /// return content merged over all chain elements
    const ACCounter & counterWeighted() const { return mergedCounterWeighted_; }

protected:
    /// branch content: global event counter stored by EventLogger
    ACCounter * counter_;
    /// branch content: global weighted event counter stored by EventLogger
    ACCounter * counterWeighted_;

    /// merged content: global event counter stored by EventLogger
    ACCounter mergedCounter_;
    /// merged content: global weighted event counter stored by EventLogger
    ACCounter mergedCounterWeighted_;
    
    /// set counter object
    ACCounter ** linkCounter();
    /// set counter object
    ACCounter ** linkCounterWeighted();
    /// merge the contents of all 'events' (i.e. one entry for each file)
    void merge();

    /// grant access to ACTreeReader
    friend class ACTreeReader;
};

#endif
