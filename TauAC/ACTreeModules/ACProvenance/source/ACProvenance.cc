#include "../interface/ACProvenance.h"

ACProvenance::ACProvenance() {
    /// pointers to branch content have to be set to 0 to prevent from root crash
    counter_ = 0;
    counterWeighted_ = 0;
}

ACCounter ** ACProvenance::linkCounter() {
    return &counter_;
}

ACCounter ** ACProvenance::linkCounterWeighted() {
    return &counterWeighted_;
}

void ACProvenance::merge() {
    //test whether the temp objects are properly connected to a branch and filled, then merge them
    if (counter_) mergedCounter_ += *counter_;
    if (counterWeighted_) mergedCounterWeighted_ += *counterWeighted_;
}