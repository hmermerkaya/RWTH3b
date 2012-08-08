// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      ThreeProngTauCreator
// 
/**
 This class creates a kinemtic tau from the 3prong decay suggestion.
 Part of the KinematicTau package.

 @author Lars Perchalla, Philip Sauerland
 @date 2009
 */

#ifndef ThreeProngTauCreator_h
#define ThreeProngTauCreator_h


#include "RecoTauTag/KinematicTau/interface/KinematicTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

#include <TLorentzVector.h>
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>
//kinematic fit:
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include <RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h>
#include <RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h>
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
//own KinematicFit classes
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
//#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackVertexLinkKinematicConstraint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class ThreeProngTauCreator : public KinematicTauCreator
{
public:
  explicit ThreeProngTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder):KinematicTauCreator(transTrackBuilder){}
  explicit ThreeProngTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const edm::ParameterSet& cfg):KinematicTauCreator(transTrackBuilder, cfg){}
  
  // ndf depends on specific decay.
  virtual int ndf() const;
  
private:
  virtual int create(SelectedKinematicDecay &KFTau);
  bool createStartScenario(SelectedKinematicDecay &KFTau, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos);

  bool kinematicRefit(std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex & primaryVertex);
  std::pair<double,double> getTauMomentumMagnitudes(double ma1,double pa1,double M,double theta);
  RefCountedKinematicParticle unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx,std::vector<TLorentzVector> &NuGuessLV);
  RefCountedKinematicParticle virtualKinematicParticle(const TransientVertex & vtxGuess, const TLorentzVector & nuGuess);
	
};

#endif
