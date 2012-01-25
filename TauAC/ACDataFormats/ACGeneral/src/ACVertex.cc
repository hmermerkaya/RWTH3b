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
double ACVertex::vtxDistanceSignificance(const ACVertex & vtx) const {
	//distance between both vertices is calculated in units of their projected errorsum
	TMatrixDSym matrix = error() + vtx.error();
	TVector3 distance = vtx.position() - position();
	double error = projectedError(distance, matrix);
	double significance = -1.;
	if(error!=0.) significance = distance.Mag()/error;
	
	return significance;
}
double ACVertex::projectedError(const TVector3 & axis, const TMatrixDSym & error) const {
	//projects the error in direction of axis
	if(axis.Mag()==0) return 0.0;
	TVector3 unit = axis.Unit();//normalize the vect
	TVectorD dist;
	dist.ResizeTo(3);
	for(unsigned int i=0; i!=3; i++) dist[i]= unit(i);
	double similarity = error.Similarity(dist);
	if(similarity<0.) similarity = 0.;//catch NaN
	
	return sqrt(similarity);
}

// some shortcuts:
const double ACVertex::vx() const { return position().X(); }
const double ACVertex::vy() const { return position().Y(); }
const double ACVertex::vz() const { return position().Z(); }
double ACVertex::svx() const { if (error()(0,0)<.0) return -1.0; return sqrt(error()(0,0)); }
double ACVertex::svy() const { if (error()(1,1)<.0) return -1.0; return sqrt(error()(1,1)); }
double ACVertex::svz() const { if (error()(2,2)<.0) return -1.0; return sqrt(error()(2,2)); }
