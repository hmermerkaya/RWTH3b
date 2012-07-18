#ifndef ACBeamSpot_H
#define ACBeamSpot_H
/**
 Base class to store the beam spot (based on the corresponding CMSSW class).

 @author Lars Perchalla
 @date 2012
 */

#include <TVector3.h>
#include <TMatrixDSym.h>

#include "ACCommon.h"

class ACBeamSpot : public ACCommon {
public:
    virtual std::string classname() const { return "ACBeamSpot"; }
    
    ACBeamSpot();
    virtual ~ACBeamSpot() {};

    const TVector3 & position() const;
    const TMatrixDSym & error() const;
    
    /// x position
    double x0() const { return position().X(); }
    /// y position 
    double y0() const { return position().Y(); }
    /// z position 
    double z0() const { return position().Z(); }

    /// x coordinate of the beeam spot position at a given z value (it takes into account the dxdz slope)
    double x(const double z) const { return x0() + dxdz() * (z - z0()); }
    /// y coordinate of the beeam spot position at a given z value (it takes into account the dydz slope)
    double y(const double z) const { return y0() + dydz() * (z - z0()); }
    /// position of the beam spot at a given z value (it takes into account the dxdz and dydz slopes)
    const TVector3 position(const double z) const;

    double sigmaZ() const { return sigmaZ_; }
    double dxdz() const { return dxdz_; }
    double dydz() const { return dydz_; }
    double beamWidthX() const { return beamWidthX_; }
    double beamWidthY() const { return beamWidthY_; }
    /// error on x
    double x0Error() const { return sqrt( error_(0,0) ); }
    /// error on y
    double y0Error() const { return sqrt( error_(1,1) ); }
    /// error on z
    double z0Error() const { return sqrt( error_(2,2) ); }
    /// error on sigma z
    double sigmaZ0Error() const { return sqrt ( error_(3,3) ); }
    /// error on dxdz
    double dxdzError() const { return sqrt ( error_(4,4) ); }
    /// error on dydz
    double dydzError() const { return sqrt ( error_(5,5) ); }    
    /// error on beam width X, assume error in X = Y
    double beamWidthXError() const { return sqrt ( error_(6,6) );}
	/// error on beam width Y, assume error in X = Y
    double beamWidthYError() const { return sqrt ( error_(6,6) );}
    int type() const { return type_; }
    
    double emittanceX() const { return emittanceX_; }
	double emittanceY() const { return emittanceY_; }
	double betaStar() const { return betaStar_; }
	double beamWidthFromBeta( double z, double e ) const;
    
    /// print common beam spot parameters
    void print( std::stringstream& ss ) const;
    
    ClassDef(ACBeamSpot,1);

protected:
    /// spatial position
    TVector3 position_;
    /// error matrix 7x7 (x, y, z, sigmaZ, dxdz, dydz, beamWidthX), assuming sigma_beamWidthX = sigma_beamWidthY
    TMatrixDSym error_;
    /// sigma z
    double sigmaZ_;
    /// beam width in x
    double beamWidthX_;
    /// beam width in y
    double beamWidthY_;
    /// dxdz slope
    double dxdz_;
    /// dydz slope
    double dydz_;
    /// emittance in x
    double emittanceX_;
    /// emittance in y
    double emittanceY_;
    /// beta star
    double betaStar_;
    /// type of current beam: Unknown=-1, Fake=0, LHC=1, Tracker=2
    int type_;
};

typedef ACBeamSpot* pACBeamSpot;

#endif
