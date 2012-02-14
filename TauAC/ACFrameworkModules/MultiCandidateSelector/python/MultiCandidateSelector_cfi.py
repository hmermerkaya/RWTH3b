import FWCore.ParameterSet.Config as cms

MultiCandidateSelector = cms.EDProducer('MultiCandidateSelector',
    filters = cms.VPSet(
        cms.PSet(
            type = cms.string("pfJets"),
            src = cms.InputTag("ak5PFJets"),
            cut = cms.string("pt > 10")
        ),
#        cms.PSet(
#            type = cms.string("muons"),
#            src = cms.InputTag("muons"),
#            cut = cms.string("pt > 10")
#        ),
#        cms.PSet(
#            type = cms.string("electrons"),
#            src = cms.InputTag("gsfElectrons"),
#            cut = cms.string("pt > 10")
#        ),
    )
)
