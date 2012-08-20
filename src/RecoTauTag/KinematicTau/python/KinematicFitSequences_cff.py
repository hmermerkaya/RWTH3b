import FWCore.ParameterSet.Config as cms
from RecoTauTag.KinematicTau.InputTrackSelector_cfi import *
from RecoTauTag.KinematicTau.ThreeProngInputSelector_Step2_cfi import *
from RecoTauTag.KinematicTau.kinematictau_cfi import *
from RecoTauTag.KinematicTau.KinematicTauSkim_cfi import *
from RecoTauTag.KinematicTau.KinematicTauAnalyzer_cfi import *
from RecoTauTag.KinematicTau.KinematicTauAnalyzer_cfi import *
from RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi import *
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *

nTauPerVtx = 1

ListofVertices = []
reducedPrimaryVtx1 = offlinePrimaryVertices
ListofVertices.append(reducedPrimaryVtx1)
if nTauPerVtx == 1:
    reducedPrimaryVtx2 = offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx2)
    reducedPrimaryVtx3 = offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx3)
    reducedPrimaryVtx4 = offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx4)
    reducedPrimaryVtx5= offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx5)
    reducedPrimaryVtx6= offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx6)
    reducedPrimaryVtx7= offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx7)
    reducedPrimaryVtx8= offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx8)
    reducedPrimaryVtx9= offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx9)
    reducedPrimaryVtx10 = offlinePrimaryVertices.clone()
    ListofVertices.append(reducedPrimaryVtx10)
    

VertexSequence   = cms.Sequence()
VertexTags       = cms.untracked.vstring()
NonTauTracksList = cms.untracked.vstring()
VertexSequences  = [cms.Sequence()]

index=0
for item in ListofVertices:
    index+=1
    NonTauTracks="NonTauTracks"
    NonTauTracks+=str(index) 
    NonTauTracksList.append(NonTauTracks)
    item.TrackLabel = cms.InputTag("InputTrackSelector",NonTauTracks)
    name+=str(index)
    print name
    VertexTags.append(name)
    currentSequence = cms.Sequence(VertexSequences[0]*item)
    VertexSequences[0] = currentSequence

VertexSequence=VertexSequences[0]

InputTrackSelector.NonTauTracks = NonTauTracksList
InputTrackSelector.nTauPerVtx = cms.untracked.uint32(nTauPerVtx)

ThreeProngInputSelectorStep2.VertexTags = VertexTags
ThreeProngInputSelectorStep2.NonTauTracks = NonTauTracksList
    
#define sequences for Kinematic Fit with Classic single vertex
KinematicFitSequence         = cms.Sequence(InputTrackSelector*VertexSequence*ThreeProngInputSelectorStep2*KinematicTauProducer)
KinematicFitSequencewithSkim = cms.Sequence(InputTrackSelector*VertexSequence*ThreeProngInputSelectorStep2*KinematicTauProducer*KinematicTauSkim)
KinematicFitSequencewithDQM  = cms.Sequence(InputTrackSelector*VertexSequence*ThreeProngInputSelectorStep2*KinematicTauProducer*KinematicTauAnalyzer)
