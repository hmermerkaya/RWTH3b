import FWCore.ParameterSet.Config as cms

process = cms.Process("InputTrackSelectorVBFH")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("HiggsKinTau.FnlAnlzr.MessageLogger_cfi")

###############
#lp parameters#
numberOfEvents = 10
verbosity = 2
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

process.load("HiggsKinTau.InputTrackSelector.InputTrackSelector_cfi")
process.InputTrackSelectorVBFH.verbosity = cms.untracked.int32(verbosity)

#process.p = cms.Path(process.tauSelectorSeq)
process.p = cms.Path(process.InputTrackSelectorVBFH)
