import FWCore.ParameterSet.Config as cms

process = cms.Process("InputTrackSelector")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_3XY_V18::All'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('KinematicTauCreator')
process.MessageLogger.debugModules = cms.untracked.vstring('InputTrackSelector')
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	KinematicTauCreator = cms.untracked.PSet(limit = cms.untracked.int32(-1))
)

###############
numberOfEvents = 10

###############
if numberOfEvents == -1:
	numberOfEvents = 10000
inPath = '/disk1/perchalla/data/CMSSW_3_1_1/KinTau/tau3piFromZ/'
jobName = 'AODSIMHLT_tau3piFromZ_10000evts'

#print 'file://'+inPath+jobName+'.root';

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file://'+inPath+jobName+'.root'),
	noEventSort = cms.untracked.bool(True), #Events are processed in order of run number and luminosity block number always. If this parameter is true, then within a luminosity block the events will be processed in the order they appear in the ROOT file TTree. If this parameter is false, then within a luminosity block the events will be processed in event number order. Defaults to false.
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")

#process.p = cms.Path(process.tauSelectorSeq)
process.p = cms.Path(process.InputTrackSelector)
