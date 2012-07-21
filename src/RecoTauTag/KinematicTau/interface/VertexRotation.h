// -*- C++ -*-
//
// Package:    VertexRotation
// Class:      VertexRotation
// 
/**
 This class should only be called from within the KinematicTau package.
 It tries to rotate the vertex link to fit into the kinematically allowed bounds of the GJ angle.
 Part of the KinematicTau package.

 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#ifndef VertexRotation_h
#define VertexRotation_h

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

	/**
	 rotates both primary and secondary vertex
	 */
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
            edm::LogError("VertexRotation")<<"ERROR! Rotation failed! "<<test;
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
		if(verbosity_>=2) LogTrace("VertexRotation")<<"movement [sigma] = "<<movement_;
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
			if(verbosity_>=1) LogTrace("VertexRotation")<<"Correction found. reached thetaMax ("<<thetaMax<<") up to "<<thetaMax-theta<<" rad.";
			success_ = true;
		}else{
			if(verbosity_>=1) LogTrace("VertexRotation")<<"Correction not found as "<<sqrt(pow(projectedErrorP,2.0)+pow(projectedErrorS,2.0))<<" < "<<correction.Mag()<<". Only reached thetaMax ("<<thetaMax<<") up to "<<thetaMax-theta<<" rad.";
			success_ = false;
		}
		
		if(verbosity_>=1){
			TVector3 pvE(pVtx.xError(), pVtx.yError(), pVtx.zError());
			TVector3 svE(sqrt(sVtx.positionError().cxx()), sqrt(sVtx.positionError().cyy()), sqrt(sVtx.positionError().czz()));
			LogTrace("VertexRotation")<<"Moved pv from ("<<pv.X()<<","<<pv.Y()<<","<<pv.Z()<<")pm("<<pvE.X()<<", "<<pvE.Y()<<", "<<pvE.Z()<<") to ("<<npv.X()<<", "<<npv.Y()<<", "<<npv.Z()<<")";
			LogTrace("VertexRotation")<<"Moved sv from ("<<sv.X()<<","<<sv.Y()<<","<<sv.Z()<<")pm("<<svE.X()<<", "<<svE.Y()<<", "<<svE.Z()<<") to ("<<nsv.X()<<", "<<nsv.Y()<<", "<<nsv.Z()<<")";
		}
		
		pVtx = newPrimVertex(npv, pVtx);
		sVtx = newSecVertex(nsv, sVtx);
		
		valid_ = true;
		
		return success_;
	}
	/**
	 rotates primary vertex around secondary one and returns the vertex significance of the modification (PV-PVrot [sigma])
	 */
	double rotatePV(reco::Vertex & pVtx, const TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir){
		double significance = 0.0;
		
		TVector3 pv(pVtx.x(), pVtx.y(), pVtx.z());
		TVector3 sv(sVtx.position().x(), sVtx.position().y(), sVtx.position().z());
		
		TVector3 ps = sv - pv;
		theta = unsignedAngle(ps, a1_.Vect());
		double thetaMax = fabs(calcThetaMax());//can only be negative for too heavy a1
		
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
	reco::Vertex newPrimVertex(TVector3 & newPoint, reco::Vertex & oldVtx){
		const math::XYZPoint point(newPoint.X(), newPoint.Y(), newPoint.Z());
		reco::Vertex vtx(point, oldVtx.error());
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
		double argument = (-pow(ma1,2.) + pow(Mtau,2.))/(2.*Mtau*pa1);// can be negative
		//catch nan
		if(fabs(argument) >  1.0) LogTrace("VertexRotation")<<"calcThetaMax: Warning! arcsin("<<argument<<") = "<<asin(argument)<<". (pa1 "<<pa1<<", ma1 "<<ma1<<")";
		if(argument >  1.0) argument =  1.0;
		if(argument < -1.0) argument = -1.0;
		
		return asin(argument);
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
		movement_ = 0.0;
	}

	double vtxDistanceSignificance(const TVector3 & pv, const TMatrixDSym & pvE, const TVector3 & sv, const TMatrixDSym & svE, TVector3 * correction = 0){
		//if correction=0 the distance between both vertices is calculated in units of their projected errorsum
		//otherwise it expresses correction in units of the projected errorsum of the vertices
		TMatrixDSym matrix = pvE + svE;
		TVector3 corr;
		if(correction==0) corr = sv - pv;
		else corr = *correction;
		double error = projectedError(corr, matrix);
		double significance = -1.;
		if(error!=0.) significance = corr.Mag()/error;
		
		return significance;
	}
	double projectedError(const TVector3 & axis, const TMatrixDSym & error){
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
	
	
};
#endif
