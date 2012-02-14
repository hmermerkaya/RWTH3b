#ifndef ACVertex_H
#define ACVertex_H
/**
 Base class to store vertex information.

 @author Lars Perchalla & Philip Sauerland
 @date 2010
 */

#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>

#include "ACCommon.h"
#include "ACRef.h"

class ACVertex : public ACCommon {
public:
    virtual std::string classname() const { return "ACVertex"; }
    ACVertex();
    ACVertex(const TVector3 & position, bool valid, const TMatrixDSym & error = TMatrixDSym(3), float chi2 = -1000., float ndof = -1000., const double & trackSize = -1000., const double & ptSum = -1000.);
    virtual ~ACVertex() {};

    const TVector3 & position() const;
    bool valid() const;
    const TMatrixDSym & error() const;
    float chi2() const;
    float ndof() const;
    double trackSize() const;
    double ptSum() const;
    /// distance from vtx in units of their projected errorsum
    double vtxDistanceSignificance(const ACVertex & vtx) const;
	/// project error in direction of axis
    double projectedError(const TVector3 & axis, const TMatrixDSym & error) const;
    
    /// position component
    const double vx() const;
    /// position component
    const double vy() const;
    /// position component
    const double vz() const;
    /// corresponding diagonal element of the covariance matrix
    double svx() const;
    /// corresponding diagonal element of the covariance matrix
    double svy() const;
    /// corresponding diagonal element of the covariance matrix
    double svz() const;
    
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
    /// ptsum of all assigned tracks (according to the official PV sorting at VertexHigherPtSquared, ignore tracks with errors larger than pt)
    double ptSum_;
};

typedef ACVertex* pACVertex;

#endif
