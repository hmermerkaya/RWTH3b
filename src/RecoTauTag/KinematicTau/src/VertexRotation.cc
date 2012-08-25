#include  "RecoTauTag/KinematicTau/interface/VertexRotation.h"

VertexRotation::VertexRotation(TLorentzVector & a1, int verbosity):
  a1_(a1), 
  verbosity_(verbosity)
{
  clear();
}

VertexRotation::VertexRotation(int verbosity):
  a1_(0,0,0,0)
{
  clear();
}

  
VertexRotation::~VertexRotation(){
}

/**
   rotates primary vertex around secondary one and returns the vertex significance of the modification (PV-PVrot [sigma])
*/
double VertexRotation::rotatePV(reco::Vertex & pVtx, const TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir,
				double reduceThetamax){
  if(a1_==TLorentzVector(0,0,0,0))return 0.0;
  double significance = 0.0;
  
  TVector3 pv(pVtx.x(), pVtx.y(), pVtx.z());
  TVector3 sv(sVtx.position().x(), sVtx.position().y(), sVtx.position().z());
  
  TVector3 ps = sv - pv;
  theta = unsignedAngle(ps, a1_.Vect());
  double thetaMax = fabs(calcThetaMax())*reduceThetamax;//can only be negative for too heavy a1
  TVector3 norm = ps.Cross(a1_.Vect());//norm vect of surface in which theta is defined
  if(norm.Mag() == 0.0) return significance;
  TRotation rot;
  rot.Rotate(theta - thetaMax, norm);//rotate right if theta > thetaMax, if not rotate left
  TVector3 psRot = rot * ps;
  TVector3 pvRot = sv - psRot;//new primary vertex position
  //now project the errors into direction of correction		
  TMatrixDSym matrixP(3);
  matrixP.ResizeTo(TMatrixDSym(3));
  for(int i=0; i!=3; i++)	for(int j=0; j!=3; j++) matrixP(i,j) = pVtx.covariance(i,j);//diagonals are squares of sigmas
  significance = vtxDistanceSignificance(pv, matrixP, pvRot, matrixP);//modification in units of sigma
  //		if(verbosity_>=2) LogTrace("VertexRotation")<<"significance [sigma] = "<<significance;
  tauFlghtDir = psRot;
  theta = unsignedAngle(tauFlghtDir, a1_.Vect());
  pVtx = newPrimVertex(pvRot, pVtx);
  return significance;
}

/**
   returns fake vertex at modified position
*/
reco::Vertex VertexRotation::newPrimVertex(TVector3 & newPoint, reco::Vertex & oldVtx){
  const math::XYZPoint point(newPoint.X(), newPoint.Y(), newPoint.Z());
  reco::Vertex vtx(point, oldVtx.error());
  return vtx;
}

TransientVertex VertexRotation::newSecVertex(TVector3 & newPoint, TransientVertex & oldVtx){
  GlobalPoint pos(newPoint.X(), newPoint.Y(), newPoint.Z());
  GlobalError posError(oldVtx.positionError());
  std::vector<reco::TransientTrack> tracks = oldVtx.originalTracks();//refittedTracks()???
  TransientVertex vtx(pos, posError, tracks, oldVtx.totalChiSquared());
  return vtx;
}

double VertexRotation::movement(){return movement_;}

bool VertexRotation::isValid(){return valid_;}

double VertexRotation::unsignedAngle(const TVector3& v1, const TVector3& v2){//smallest angle between v1 and v2. uses TVector3.Angle(v2)
  double angle = v1.Angle(v2);
  while(angle >= TMath::Pi())       angle -= (2.0*TMath::Pi());
  while(angle < (-1.0*TMath::Pi())) angle += (2.0*TMath::Pi());
  return fabs(angle);
}

double VertexRotation::calcThetaMax(){
  double ma1 = a1_.M(), pa1 = a1_.P(), Mtau = 1.777;
  double argument = (-pow(ma1,2.) + pow(Mtau,2.))/(2.*Mtau*pa1);// can be negative
  //catch nan
  if(fabs(argument) >  1.0) LogTrace("VertexRotation")<<"calcThetaMax: Warning! arcsin("<<argument<<") = "<<asin(argument)<<". (pa1 "<<pa1<<", ma1 "<<ma1<<")";
  if(argument >  1.0) argument =  1.0;
  if(argument < -1.0) argument = -1.0;
  return asin(argument);
}
  
	
void VertexRotation::clear(){
  valid_ = false;
  success_ = false;
  movement_ = 0.0;
}
  
double VertexRotation::vtxDistanceSignificance(TVector3 & pv,TMatrixDSym & pvE,TVector3 & sv, TMatrixDSym & svE){
  //if correction=0 the distance between both vertices is calculated in units of their projected errorsum
  //otherwise it expresses correction in units of the projected errorsum of the vertices
  TMatrixDSym matrix = pvE + svE;
  TVector3 corr = sv - pv;
    double error = projectedError(corr, matrix);
  double significance = -1.;
  if(error!=0.) significance = corr.Mag()/error;
  return significance;
}

double VertexRotation::projectedError(const TVector3 & axis, const TMatrixDSym & error){
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
	
