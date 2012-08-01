import FWCore.ParameterSet.Config as cms
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *

reducedPrimaryVertices = offlinePrimaryVertices
reducedPrimaryVertices.TrackLabel = cms.InputTag("InputTrackSelector", "NonTauTracks")

