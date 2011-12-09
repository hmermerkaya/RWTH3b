#include "../interface/ACVertex.h"

ClassImp(ACVertex)

ACVertex::ACVertex() {
    position_.SetXYZ(-1000., -1000., -1000.);
    error_.ResizeTo(TMatrixDSym(3));
    for (int i=0; i!=3; i++)    for (int j=0; j!=3; j++) error_(i,j) = -1000.0;
    valid_ = false;
    chi2_ = -1000.;
    ndof_ = -1000.;
    trackSize_ = -1000;
//    ptSum_ = -1000.;
}
ACVertex::ACVertex(const TVector3 & position, bool valid, const TMatrixDSym & error, float chi2, float ndof, int trackSize)://double ptSum
position_(position), valid_(valid), chi2_(chi2), ndof_(ndof), trackSize_(trackSize)
// , ptSum_(ptSum)
{
    error_.ResizeTo(TMatrixDSym(3));
    error_ = error;
}
TVector3 ACVertex::position() const { return position_; }
bool ACVertex::valid() const { return valid_; }
TMatrixDSym ACVertex::error() const { return error_; }
float ACVertex::chi2() const { return chi2_; }
float ACVertex::ndof() const { return ndof_; }
int ACVertex::trackSize() const { return trackSize_; }
// double ptSum() const { return ptSum_; }
