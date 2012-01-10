import FWCore.ParameterSet.Config as cms
import os
import ACFrameworkModules.Common.Tools as Tools
import ACFrameworkModules.FinalTreeFiller.Steering as Steering

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

process = cms.Process("FinalTreeFiller")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("ACFrameworkModules.Common.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


###############
#steering parameters#
inputPath, outputPath, jobName, isData, jsonFile = Steering.LoadUserParameters()
numberOfEvents = -1
ignoreFilter = True # but always ignore TauMotherProducer until its ready
printEvents = 0	#print generator event
triggerTag = "HLT"
minTau = 1 #minimum of selected taus
###############

### scan all files located in inputPath
files = Tools.scanDir(inputPath)
sourceString = Tools.buildPoolSource(files)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(sourceString)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

### select runs from JSON file
if isData:
    print 'Dicriminate lumisections from', jsonFile
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()#convert into compact format needed by CMSSW
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(lumiString.split(','))

if isData:
    process.GlobalTag.globaltag = 'GR_P_V27::All'
else:
    process.GlobalTag.globaltag = 'START44_V7::All'

print 'Using the GlobalTag ',process.GlobalTag.globaltag

motherPdgID = 0
if isData:
    decayType = 'data'
else:
    decayType, motherPdgID = Tools.parseDecayType(jobName)

print 'Analysing',numberOfEvents,'events in <'+inputPath+'>. Detected',len(sourceString) ,'files of decayType = '+decayType+'(',motherPdgID,').'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('file://'+outputPath+jobName+'_'+str(numberOfEvents)+'_7TeV.root'), #+date.today().strftime('%Y%m%d')+'.root'),
    closeFileFast = cms.untracked.bool(True)
)

### generator workflow
process.load("ACFrameworkModules.GenSelector.GenSelector_cfi")
process.printTree.maxEventsToPrint = cms.untracked.int32(printEvents)
process.printSelectedTree.maxEventsToPrint = cms.untracked.int32(printEvents)
process.GenSelector.decayType = cms.untracked.string(decayType)
process.GenSelector.motherPdgID = cms.untracked.int32(motherPdgID)

### kinematic fit workflow
process.load("CommonTools.PrimVtxSelector.PrimVtxSelector_cfi")
process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")
process.InputTrackSelector.minTau = cms.untracked.uint32(minTau)
process.InputTrackSelector.filterTaus = cms.untracked.bool(False)
process.load("RecoTauTag.KinematicTau.ThreeProngInputSelector_cfi")
process.ThreeProngInputSelector.minTau = cms.untracked.uint32(minTau)
process.load("RecoTauTag.KinematicTau.kinematictauAdvanced_cfi")
process.KinematicTauProducer.minKinTau = cms.untracked.uint32(minTau)
process.KinematicTauProducer.primVtx = cms.InputTag("ThreeProngInputSelector","primVtx")
process.KinematicTauProducer.selectedTauCandidates = cms.InputTag("ThreeProngInputSelector","InputTauRefs")
process.KinematicTauProducer.inputTracks = cms.InputTag("ThreeProngInputSelector","InputTracks")

### event storage
process.load("ACFrameworkModules.FinalTreeFiller.FinalTreeFiller_cfi")
process.FinalTreeFiller.triggerResults = cms.InputTag("TriggerResults","",triggerTag)
process.FinalTreeFiller.decayType = cms.untracked.string(decayType)

#ignore filter
process.ignorePath = cms.Path(
    process.generatorSequence
    *cms.ignore(process.PrimVtxSelector)
    #*cms.ignore(process.HLTSelector)
    *cms.ignore(process.InputTrackSelector)
    *cms.ignore(process.ThreeProngInputSelector)
    *cms.ignore(process.KinematicTauProducer)
    *process.matchingSeq
    *process.FinalTreeFiller
)
#cumulative
process.cumulativePath = cms.Path(
    process.generatorSequence
    *process.PrimVtxSelector
    #*process.HLTSelector
    *process.InputTrackSelector
    *process.ThreeProngInputSelector
    *process.KinematicTauProducer
    *process.matchingSeq
    *process.FinalTreeFiller
)

#if no schedule is provided all pathes are executed
if ignoreFilter:
	process.schedule = cms.Schedule(process.ignorePath)
else:
	process.schedule = cms.Schedule(process.cumulativePath)

