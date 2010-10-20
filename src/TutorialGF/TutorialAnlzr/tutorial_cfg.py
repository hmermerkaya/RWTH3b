import FWCore.ParameterSet.Config as cms

process = cms.Process("TutorialGF")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.GlobalTag.globaltag = 'START36_V9::All'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0))
)

numberOfEvents = 100 #-1

#specify input file
inPath = '/disk1/perchalla/data/CMSSW_3_6_2/ZTauTau3pr-START36_V9/AOD-HLT8E29/0544449c08827763af515821e3b042e2/'
jobName = 'AODHLT_tau3piFromZ_7000GeV_1000evts_1'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file://'+inPath+jobName+'.root'),
	noEventSort = cms.untracked.bool(True), #Events are processed in order of run number and luminosity block number always. If this parameter is true, then within a luminosity block the events will be processed in the order they appear in the ROOT file TTree. If this parameter is false, then within a luminosity block the events will be processed in event number order. Defaults to false.
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

#load the module definition
process.load("TutorialGF.GenSelectorGF.GenSelectorGF_cfi")
process.load("TutorialGF.TutorialAnlzr.TutorialAnlzr_cfi")

#define the path for the event loop
process.p = cms.Path(process.GenSelector*process.TutorialAnlzr)
