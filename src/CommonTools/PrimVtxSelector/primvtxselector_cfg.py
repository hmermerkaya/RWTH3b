import FWCore.ParameterSet.Config as cms

process = cms.Process("PrimVtxSelectorVBFH")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.debugModules = cms.untracked.vstring('PrimVtxSelector')
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1))
)

numberOfEvents = -1

inPath = '/disk1/perchalla/data/CMSSW_3_1_2/KinTau/tau3piFromVBFH/'
jobName = 'AODSIMHLT_tau3piFromVBFH_145GeV'

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

process.load("CommonTools.PrimVtxSelector.PrimVtxSelector_cfi")

process.p = cms.Path(process.PrimVtxSelector)
