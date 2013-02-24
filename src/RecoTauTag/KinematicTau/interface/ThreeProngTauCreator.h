// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      ThreeProngTauCreator
// 
/*
 @author Lars Perchalla, Philip Sauerland
 @date 2009
 Modified by Ian M. Nugent
*/

#ifndef ThreeProngTauCreator_h
#define ThreeProngTauCreator_h


#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

#include <TLorentzVector.h>
//own KinematicFit classes
#include "RecoTauTag/KinematicTau/interface/ParticleMassHelper.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "RecoTauTag/KinematicTau/interface/FitSequencer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


class ThreeProngTauCreator : public FitSequencer{
public:
  enum FitSeq{VertexFit,TauFit,NFits};
  
  explicit ThreeProngTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder,edm::Handle<reco::GenParticleCollection> &GenPart_):FitSequencer(transTrackBuilder,GenPart_){}
  explicit ThreeProngTauCreator(edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const edm::ParameterSet& cfg,edm::Handle<reco::GenParticleCollection> &GenPart_):FitSequencer(transTrackBuilder, cfg, GenPart_){}
  
  TString FitSequence(int i){
    if(i==VertexFit) return "VertexFit";
    if(i==TauFit)    return "TauFit";
    return "Invalid";
  }

private:
  int  create(unsigned int& ambiguity,SelectedKinematicDecay &KFTau);
  void ConfigurePions(SelectedKinematicDecay &KFTau, std::vector<TrackParticle> &pions);
  bool FitA1(std::vector<TrackParticle> &pions,const reco::Vertex & primaryVertex);
  bool FitTau(std::vector<LorentzVectorParticle>  &unfitDaughters,const reco::Vertex & primaryVertex,unsigned int &ambiguity);

  // Parameters
  ParticleMassHelper PMH;

};

#endif
