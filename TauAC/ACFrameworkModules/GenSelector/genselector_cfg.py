import FWCore.ParameterSet.Config as cms

process = cms.Process("GenSelector")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("HiggsKinTau.FnlAnlzr.MessageLogger_cfi")

###############
#lp parameters#
numberOfEvents = 10
printEvents = 1
verbosity = 2
###############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)
inPath = '/disk1/perchalla/data/CMSSW_3_6_2/ZTauTau3pr-START36_V9/AOD-HLT8E29/0544449c08827763af515821e3b042e2/'
jobName = 'AODHLT_tau3piFromZ_7000GeV_1000evts_1'

#print 'file://'+inPath+jobName+'.root';

decayType = 'unknown';
if jobName.find('tau3piFromVBFH')!=-1:
	decayType = 'VBFH'
if jobName.find('tau3piFromZ')!=-1:
	decayType = 'Z3pr'
if jobName.find('DYTauTau')!=-1:
	decayType = 'ZSMBr'

if numberOfEvents == -1:
	numberOfEvents = 10000

print 'Analysing',numberOfEvents,'events in <'+inPath+jobName+'>. Detected decayType = '+decayType+'.'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file://'+inPath+jobName+'.root'),
	noEventSort = cms.untracked.bool(True), #Events are processed in order of run number and luminosity block number always. If this parameter is true, then within a luminosity block the events will be processed in the order they appear in the ROOT file TTree. If this parameter is false, then within a luminosity block the events will be processed in event number order. Defaults to false.
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.load("HiggsKinTau.GenSelector.GenSelector_cfi")
process.printTree.maxEventsToPrint = cms.untracked.int32(printEvents)
process.printSelectedTree.maxEventsToPrint = cms.untracked.int32(printEvents)
process.GenSelector.decayType = cms.untracked.string(decayType)
process.GenSelector.verbosity = cms.untracked.int32(verbosity)

process.p = cms.Path(process.genSeq*process.printTree*process.printSelectedTree)
