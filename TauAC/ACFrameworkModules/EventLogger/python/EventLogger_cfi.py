import FWCore.ParameterSet.Config as cms

EventLogger = cms.EDAnalyzer("EventLogger",
    triggerResults = cms.InputTag('TriggerResults','','FinalTreeFiller'),
    pathName = cms.string('cumulativePath'),    #pathname to be logged
    modules = cms.untracked.vstring(),        #if non empty, only modules with matching labels are logged
    #set the following for an additional set of weighted counters (MC only)
    pileUpDistributionFileMC = cms.untracked.string(""),
    pileUpDistributionHistMC = cms.untracked.string("pileup"),
    pileUpDistributionFileData = cms.untracked.string(""),
    pileUpDistributionHistData = cms.untracked.string("pileup"),
)
