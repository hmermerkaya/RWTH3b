#include "../interface/ACBeamSpot.h"

ClassImp(ACBeamSpot)

ACBeamSpot::ACBeamSpot() {
    position_.SetXYZ(-1000., -1000., -1000.);
    error_.ResizeTo(TMatrixDSym(7));
    for (int i=0; i!=7; i++) for (int j=0; j!=7; j++) error_(i,j) = -1000.0;
    sigmaZ_ = -1.;
    beamWidthX_ = -1.;
    beamWidthY_ = -1.;
    dxdz_ = -1.;
    dydz_ = -1.;
    emittanceX_ = -1.;
    emittanceY_ = -1.;
    betaStar_ = -1.;
    type_ = -1;
}

const TVector3 & ACBeamSpot::position() const { return position_; }
const TMatrixDSym & ACBeamSpot::error() const { return error_; }

const TVector3 ACBeamSpot::position(const double z) const {
    TVector3 pos(x(z), y(z), z);
    return pos;
}
double ACBeamSpot::beamWidthFromBeta( double z, double e ) const {
    return sqrt( e*betaStar_*(1 + pow((z-position_.Z())/betaStar_,2) ) );
}
void ACBeamSpot::print(std::stringstream& ss) const {
    ss << "-----------------------------------------------------\n"
    << "              Beam Spot Data\n\n"
    << " Beam type    = " << type() << "\n"
    << "       X0     = " << x0() << " +/- " << x0Error() << " [cm]\n"
    << "       Y0     = " << y0() << " +/- " << y0Error() << " [cm]\n"
    << "       Z0     = " << z0() << " +/- " << z0Error() << " [cm]\n"
    << " Sigma Z0     = " << sigmaZ() << " +/- " << sigmaZ0Error() << " [cm]\n"
    << " dxdz         = " << dxdz() << " +/- " << dxdzError() << " [radians]\n"
    << " dydz         = " << dydz() << " +/- " << dydzError() << " [radians]\n"
    << " Beam width X = " << beamWidthX() << " +/- " << beamWidthXError() << " [cm]\n"
    << " Beam width Y = " << beamWidthY() << " +/- " << beamWidthYError() << " [cm]\n"
    << " EmittanceX   = " << emittanceX() << " [cm]\n"
    << " EmittanceY   = " << emittanceY() << " [cm]\n"
    << " beta-star    = " << betaStar() << " [cm]\n"
    << "-----------------------------------------------------\n\n";
}
