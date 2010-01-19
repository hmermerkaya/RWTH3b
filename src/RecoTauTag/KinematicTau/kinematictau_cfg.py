import FWCore.ParameterSet.Config as cms

process = cms.Process("KinTauCreator")

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

process.load("HiggsKinTau.PrimVtxSelector.PrimVtxSelector_cfi")
process.PrimVtxSelectorVBFH.verbosity = cms.untracked.int32(verbosity)
process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")
process.InputTrackSelector.verbosity = cms.untracked.int32(verbosity)
process.load("RecoTauTag.KinematicTau.kinematictau_cfi")
process.KinematicTauProducer.verbosity = cms.untracked.int32(verbosity)

#process.p = cms.Path(process.tauSelectorSeq)
process.p = cms.Path(process.PrimVtxSelectorVBFH*process.InputTrackSelector*process.KinematicTauProducer)
