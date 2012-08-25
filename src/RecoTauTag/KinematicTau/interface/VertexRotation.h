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
  VertexRotation(int verbosity = 0);
  ~VertexRotation();
  
  double rotatePV(reco::Vertex & pVtx, const TransientVertex & sVtx, double & theta, TVector3 & tauFlghtDir,double reduceThetamax=1.0);
  reco::Vertex newPrimVertex(TVector3 & newPoint, reco::Vertex & oldVtx);  
  TransientVertex newSecVertex(TVector3 & newPoint, TransientVertex & oldVtx);
  double movement();
  bool isValid();
  double unsignedAngle(const TVector3& v1, const TVector3& v2);
  double calcThetaMax();
  double vtxDistanceSignificance(TVector3 & pv,TMatrixDSym & pvE,TVector3 & sv, TMatrixDSym & svE);
	
private:
  bool valid_, success_;
  TLorentzVector a1_;
  double movement_;
  int verbosity_;
  
  void clear();  
  double projectedError(const TVector3 & axis, const TMatrixDSym & error);	
	
};
#endif
