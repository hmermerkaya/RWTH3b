#ifndef ACVertex_H
#define ACVertex_H
/**
 Base class to store vertex information.

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */

#include <TVector3.h>
#include <TMatrixDSym.h>

#include "ACCommon.h"
#include "ACRef.h"

class ACVertex : public ACCommon {
public:
    virtual std::string classname() const { return "ACVertex"; }
    ACVertex();
    ACVertex(const TVector3 & position, bool valid, const TMatrixDSym & error = TMatrixDSym(3), float chi2 = -1000., float ndof = -1000., int trackSize = -1000); //double ptSum = -1000.
    virtual ~ACVertex() {};

    TVector3 position() const;
    bool valid() const;
    TMatrixDSym error() const;
    float chi2() const;
    float ndof() const;
    int trackSize() const;
//    double ptSum() const;

    ClassDef(ACVertex,1);

protected:
    /// spatial position
    TVector3 position_;
    /// whether the vertex is valid
    bool valid_;
    /// error matrix
    TMatrixDSym error_;
    /// chi2 of the vertex fit
    float chi2_;
    /// ndof of the vertex fit
    float ndof_;
    /// size of all assigned tracks
    double trackSize_;
//    /// ptsum of all assigned tracks
//    double ptSum_;
};

typedef ACVertex* pACVertex;

#endif
