// -*- C++ -*-
//
// Package:    VertexRotation
// Class:      VertexRotation
// 
/**
 
 Description: tries to rotate vertex link to fit into kinematically allowd bounds of the GJ angle
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Thu Dec  16 11:12:54 CEST 2009
// $Id: VertexRotation.h,v 1.3 2010/01/27 14:28:53 perchall Exp $
//
//

// system include files
#include <memory>
#include <fstream>

#include <TMatrixDSym.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
//vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>


class VertexRotation{
public:
	VertexRotation(TLorentzVector & a1, int verbosity = 0):
	a1_(a1), verbosity_(verbosity)
	{
		clear();
	};
	~VertexRotation(){
	};


	bool tryCorrection(reco::Vertex & pVtx, TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir, bool forceRotation = true){
		valid_ = false, success_ = false;
		TVector3 pv(pVtx.x(), pVtx.y(), pVtx.z());
		TVector3 sv(sVtx.position().x(), sVtx.position().y(), sVtx.position().z());
		
		TVector3 ps = sv - pv;
		theta = unsignedAngle(ps, a1_.Vect());
		double thetaMax = fabs(calcThetaMax());//can only be negative for too heavy a1
		
		TVector3 norm = ps.Cross(a1_.Vect());//norm vect of surface in which theta is defined
		if(norm.Mag() == 0.0) return false;
		TRotation rot;
		rot.Rotate(theta - thetaMax, norm);//rotate right if theta > thetaMax, if not rotate left
		TVector3 psRot = rot * ps;
		double test = fabs(unsignedAngle(psRot, a1_.Vect()) - thetaMax);
		if(test > pow(10., -6.)){
			printf("VertexRotation::ERROR! Rotation failed! %f\n", test);
			valid_ = false;
			return false;
		}
		TVector3 correction = ps - psRot;

		//now project the errors into direction of correction		
		TMatrixDSym matrixP(3);
		matrixP.ResizeTo(TMatrixDSym(3));
		for(int i=0; i!=3; i++)	for(int j=0; j!=3; j++) matrixP(i,j) = pVtx.covariance(i,j);//diagonals are squares of sigmas
		TMatrixDSym matrixS(3);
		matrixS.ResizeTo(TMatrixDSym(3));
		for(int i=0; i!=3; i++)	for(int j=0; j!=3; j++) matrixS(i,j) = sVtx.positionError().matrix()(i+1,j+1);//diagonals are squares of sigmas
		movement_ = vtxDistanceSignificance(pv, matrixP, sv, matrixS, &correction);//correction in units of sigma
		if(verbosity_>=2) std::cout<<"movement [sigma] = "<<movement_<<std::endl;
		TVector3 npv, nsv;
		double projectedErrorP = projectedError(correction, matrixP);
		double projectedErrorS = projectedError(correction, matrixS);
		
		if(movement_ <= 1 || forceRotation){//correction is smaller than 1sigma or rotation is forced not regarding the errors
			//allocate the full correction on both vertices according to their errors projected in correction direction
			//no need of quadratic addition?!?
//			npv = pv + (projectedErrorP/(projectedErrorP+projectedErrorS))*correction;
//			nsv = sv - (projectedErrorS/(projectedErrorP+projectedErrorS))*correction;
			npv = pv + correction;
			nsv = sv;
		}else{
			//allocate only that part of the correction on both vertices that fits into their errors
			//no need of quadratic addition?!?
			npv = pv + (projectedErrorP/correction.Mag())*correction;
			nsv = sv - (projectedErrorS/correction.Mag())*correction;
		}
		tauFlghtDir = nsv - npv;
		theta = unsignedAngle(tauFlghtDir, a1_.Vect());
		if(movement_ <= 1){
			if(verbosity_>=1) printf("VertexRotation:: correction found. reached thetaMax (%f) up to %f rad.\n", thetaMax, thetaMax-theta);
			success_ = true;
		}else{
			if(verbosity_>=1) printf("VertexRotation:: correction not found as %f < %f. Only reached thetaMax (%f) up to %f rad.\n", sqrt(pow(projectedErrorP,2.0)+pow(projectedErrorS,2.0)), correction.Mag(), thetaMax, thetaMax-theta);
			success_ = false;
		}
		
		if(verbosity_>=1){
			TVector3 pvE(pVtx.xError(), pVtx.yError(), pVtx.zError());
			TVector3 svE(sqrt(sVtx.positionError().cxx()), sqrt(sVtx.positionError().cyy()), sqrt(sVtx.positionError().czz()));
			printf("VertexRotation:: Moved pv from (%f,%f,%f)pm(%f, %f, %f) to (%f, %f, %f)\n",
				   pv.X(), pv.Y(), pv.Z(), pvE.X(), pvE.Y(), pvE.Z(), npv.X(), npv.Y(), npv.Z()
				   );
			printf("VertexRotation:: Moved sv from (%f,%f,%f)pm(%f, %f, %f) to (%f, %f, %f)\n",
				   sv.X(), sv.Y(), sv.Z(), svE.X(), svE.Y(), svE.Z(), nsv.X(), nsv.Y(), nsv.Z()
				   );
		}
		
		pVtx = newPrimVertex(npv, pVtx);
		sVtx = newSecVertex(nsv, sVtx);
		
		valid_ = true;
		
		return success_;
	}

	reco::Vertex newPrimVertex(TVector3 & newPoint, reco::Vertex & oldVtx){
		const math::XYZPoint point(newPoint.X(), newPoint.Y(), newPoint.Z());
		reco::Vertex vtx(point, oldVtx.error(), oldVtx.chi2(), oldVtx.ndof(), oldVtx.tracksSize());
		return vtx;
	}
	TransientVertex newSecVertex(TVector3 & newPoint, TransientVertex & oldVtx){
		GlobalPoint pos(newPoint.X(), newPoint.Y(), newPoint.Z());
		GlobalError posError(oldVtx.positionError());
		std::vector<reco::TransientTrack> tracks = oldVtx.originalTracks();//refittedTracks()???
		TransientVertex vtx(pos, posError, tracks, oldVtx.totalChiSquared());
		return vtx;
	}
	double movement(){return movement_;}
	bool isValid(){return valid_;}
	
	double unsignedAngle(const TVector3& v1, const TVector3& v2){//smallest angle between v1 and v2. uses TVector3.Angle(v2)
		double angle = v1.Angle(v2);
		while(angle >= TMath::Pi())       angle -= (2.0*TMath::Pi());
		while(angle < (-1.0*TMath::Pi())) angle += (2.0*TMath::Pi());
		
		return fabs(angle);
	}
	double calcThetaMax(){
		double ma1 = a1_.M(), pa1 = a1_.P(), Mtau = 1.777;
		double thetaMax = asin((-pow(ma1,2.) + pow(Mtau,2.))/(2.*Mtau*pa1));// can be negative
		
		return thetaMax;
	}
	
	/**
	 dump to ascii to be used in mathematica
	 */
	void dumpEvt(const TVector3 & pv, const TVector3 & sv, double thetaMax){
		TVector3 pvE(0,0,0);
		TVector3 svE(0,0,0);
		dumpEvt(pv, pvE, sv, svE, thetaMax);
	}
	void dumpEvt(const reco::Vertex & pVtx, const TransientVertex & sVtx, double thetaMax){
		TVector3 pv(pVtx.x(), pVtx.y(), pVtx.z());
		TVector3 sv(sVtx.position().x(), sVtx.position().y(), sVtx.position().z());
		TVector3 pvE(pVtx.xError(), pVtx.yError(), pVtx.zError());
		TVector3 svE(sVtx.positionError().cxx(), sVtx.positionError().cyy(), sVtx.positionError().czz());
		dumpEvt(pv, pvE, sv, svE, thetaMax);
	}
	void dumpEvt(const TVector3 & pv, const TVector3 & pvE, const TVector3 & sv, const TVector3 & svE, double thetaMax){
		std::ofstream fout("../../output/asciiDump/VertexRotation.dat", std::ios_base::app);//std::ios::trunc);
		//              fout.setf(std::ios::scientific);
		//              fout.setf(std::ios::fixed);//fixed
		fout<<std::setprecision(18);
		
		fout <<"<evt>\n";
		fout <<"pv";
		for(unsigned int i=0; i!=3; i++)fout<<"\t"<<pv(i);
		fout <<"\n";
		
		fout <<"pvE";
		for(unsigned int i=0; i!=3; i++)fout<<"\t"<<pvE(i);
		fout <<"\n";
		
		fout <<"sv";
		for(unsigned int i=0; i!=3; i++)fout<<"\t"<<sv(i);
		fout <<"\n";
		
		fout <<"svE";
		for(unsigned int i=0; i!=3; i++)fout<<"\t"<<svE(i);
		fout <<"\n";
		
		fout <<"a1";
		for(unsigned int i=0; i!=3; i++)fout<<"\t"<<a1_.Vect()(i);
		fout <<"\n";
		
		fout <<"thetaMax";
		fout<<"\t"<<thetaMax;
		fout <<"\n";
		fout <<"</evt>\n";
		
		fout.close();
	}
	
private:
	bool valid_, success_;
	TLorentzVector a1_;
	double movement_;
	int verbosity_;
	
	void clear(){
		valid_ = false;
		success_ = false;
	}

	double vtxDistanceSignificance(const TVector3 & pv, const TMatrixDSym & pvE, const TVector3 & sv, const TMatrixDSym & svE, TVector3 * correction = 0){
		//if correction=0 the distance between both vertices is calculated in units of their projected errorsum
		//otherwise it expresses correction in units of the projected errorsum of the vertices
		TMatrixDSym matrix = pvE + svE;
		TVector3 corr;
		if(correction==0) corr = sv - pv;
		else corr = *correction;
		double error = projectedError(corr, matrix);
		double significance = 0.;
		if(error!=0) significance = corr.Mag()/error;
		
		return significance;
	}
	double projectedError(const TVector3 & axis, const TMatrixDSym & error){
		//projects the error in direction of axis
		if(axis.Mag()==0) return 0.0;
		TVector3 unit = axis.Unit();//normalize the vect
		TVectorD dist;
		dist.ResizeTo(3);
		for(unsigned int i=0; i!=3; i++) dist[i]= unit(i);
		
		return sqrt(error.Similarity(dist));
	}
	
	
};
