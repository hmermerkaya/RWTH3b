import FWCore.ParameterSet.Config as cms
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *

reducedPrimaryVerticesNonTauTracks = offlinePrimaryVertices
reducedPrimaryVerticesNonTauTracks.TrackLabel = cms.InputTag("InputTrackSelector", "NonTauTracks")

