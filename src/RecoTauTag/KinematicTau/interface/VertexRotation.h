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
#include "FWCore/MessageLogger/interface/MessageLogger.h"


class VertexRotation {
 public:
  VertexRotation(TLorentzVector & a1, int verbosity = 0);
  ~VertexRotation();
  
  /**
     rotates both primary and secondary vertex
  */
  bool tryCorrection(reco::Vertex & pVtx, TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir, bool forceRotation = true);
  double rotatePV(reco::Vertex & pVtx, const TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir);
  reco::Vertex newPrimVertex(TVector3 & newPoint, reco::Vertex & oldVtx);  
  TransientVertex newSecVertex(TVector3 & newPoint, TransientVertex & oldVtx);
  double movement();
  bool isValid();
  double unsignedAngle(const TVector3& v1, const TVector3& v2);
  double calcThetaMax();
  void dumpEvt(const TVector3 & pv, const TVector3 & sv, double thetaMax);  
  void dumpEvt(const reco::Vertex & pVtx, const TransientVertex & sVtx, double thetaMax);
  void dumpEvt(const TVector3 & pv, const TVector3 & pvE, const TVector3 & sv, const TVector3 & svE, double thetaMax);
	
private:
  bool valid_, success_;
  TLorentzVector a1_;
  double movement_;
  int verbosity_;
  
  void clear();  
  double vtxDistanceSignificance(const TVector3 & pv, const TMatrixDSym & pvE, const TVector3 & sv, const TMatrixDSym & svE, TVector3 * correction = 0);
  double projectedError(const TVector3 & axis, const TMatrixDSym & error);	
	
};
#endif
