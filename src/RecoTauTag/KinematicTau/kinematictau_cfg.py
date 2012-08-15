import FWCore.ParameterSet.Config as cms

process = cms.Process("KinematicTauAnalyzer")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.GlobalTag.globaltag = 'START36_V9::All'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('KinematicTauCreator')
process.MessageLogger.debugModules = cms.untracked.vstring('KinematicTauProducer','InputTrackSelector')
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
	KinematicTauCreator = cms.untracked.PSet(limit = cms.untracked.int32(-1))
)

numberOfEvents = 100

inPath = '/disk1/perchalla/data/CMSSW_3_1_1/KinTau/tau3piFromZ/'
jobName = 'AODSIMHLT_tau3piFromZ_10000evts'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file://'+inPath+jobName+'.root'),
	noEventSort = cms.untracked.bool(True), #Events are processed in order of run number and luminosity block number always. If this parameter is true, then within a luminosity block the events will be processed in the order they appear in the ROOT file TTree. If this parameter is false, then within a luminosity block the events will be processed in event number order. Defaults to false.
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

process.load("RecoTauTag.KinematicTau.KinematicFitSequences_cff")

process.p = cms.Path(process.KinematicFitSequence)
