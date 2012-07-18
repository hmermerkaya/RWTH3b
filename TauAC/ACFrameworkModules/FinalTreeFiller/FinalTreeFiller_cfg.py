import FWCore.ParameterSet.Config as cms
import os
import ACFrameworkModules.Common.Tools as Tools
import ACFrameworkModules.FinalTreeFiller.Steering as Steering
import ACFrameworkModules.FinalTreeFiller.PileUpReweighting as PileUpReweighting

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

#relative path to $CMSSW_BASE/src. needed for crab and to run from arbitrary paths
baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")

processName = "FinalTreeFiller"
process = cms.Process(processName)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("ACFrameworkModules.Common.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

###############
#steering parameters#
inputPath, outputPath, jobName, globaltag, isData, jsonFile, pileUpDistributionMC = Steering.LoadUserParameters()
numberOfEvents = 100000
ignoreFilter = True # but always ignore TauMotherProducer until its ready
printEvents = 0	#print generator event
triggerTag = "HLT"
minTau = 1 #minimum of selected taus
minJetPt = 10. #ignore pfjets below this pt threshold (in GeV)
minTauPt = 10. #ignore pftaus below this pt threshold (in GeV)
minTrackPt = 1. #ignore tracks below this pt threshold (in GeV)
rerunTauID = True #whether or not to rerun the PFTau sequence
doPileUpReweighting = True
pileUpReweightingType = "true" #(true or observed) choose according to https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities
runFromCrab = True
###############

if not runFromCrab:
    ### when executed locally scan all files located in inputPath
    files = Tools.scanDir(inputPath)
    sourceString = Tools.buildPoolSource(files)
else:
    ### in grid mode define input via crab.cfg. python will be send precompiled, for the c++ part baseDir needs to be the absolute path (starting from the from the CMSSW_BASE)
    sourceString = ""
    baseDir = 'src'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(sourceString)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

dataPileUpFilename = ''
if doPileUpReweighting:
    dataPileUpFilename = baseDir + PileUpReweighting.createDataPileUpFile(jsonFile, pileUpReweightingType)

### select runs from JSON file
if isData:
    print 'Discriminate lumisections from', jsonFile
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()#convert into compact format needed by CMSSW
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(lumiString.split(','))

process.GlobalTag.globaltag = globaltag
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
if runFromCrab:
    process.TFileService.fileName = cms.string(jobName+'_'+str(numberOfEvents)+'_7TeV.root')

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
process.InputTrackSelector.minTauPt = cms.untracked.double(minTauPt)
process.load("RecoTauTag.KinematicTau.ThreeProngInputSelector_cff")
process.ThreeProngInputSelectorStep1.minTau = cms.untracked.uint32(minTau)
process.ThreeProngInputSelectorStep2.minTau = cms.untracked.uint32(minTau)
process.load("RecoTauTag.KinematicTau.kinematictauAdvanced_cfi")
process.KinematicTauProducer.minKinTau = cms.untracked.uint32(minTau)
process.KinematicTauProducer.primVtx = cms.InputTag("ThreeProngInputSelectorStep2", "primVtx")
process.KinematicTauProducer.selectedTauCandidates = cms.InputTag("ThreeProngInputSelectorStep2", "InputTauRefs")
process.KinematicTauProducer.inputTracks = cms.InputTag("ThreeProngInputSelectorStep2", "InputTracks")
process.load("RecoTauTag.KinematicTau.kinematictau_cfi")
process.load("RecoTauTag.KinematicTau.KinematicTauSkim_cfi")

### event storage
process.load("ACFrameworkModules.FinalTreeFiller.FinalTreeFiller_cfi")
process.FinalTreeFiller.triggerResults = cms.InputTag("TriggerResults","",triggerTag)
process.FinalTreeFiller.decayType = cms.untracked.string(decayType)
process.FinalTreeFiller.minTauPt = cms.untracked.double(minTauPt)
process.FinalTreeFiller.minJetPt = cms.untracked.double(minJetPt)
process.FinalTreeFiller.minTrackPt = cms.untracked.double(minTrackPt)

if doPileUpReweighting and not isData:
    process.FinalTreeFiller.pileUpDistributionFileMC = cms.untracked.string(baseDir+"/data/MCPileUpDistributions.root")
    process.FinalTreeFiller.pileUpDistributionHistMC = cms.untracked.string(pileUpDistributionMC)
    process.FinalTreeFiller.pileUpDistributionFileData = cms.untracked.string(dataPileUpFilename)
    process.FinalTreeFiller.pileUpDistributionHistData = cms.untracked.string("pileup")

### candidate selectors (deprecated as violating refs)
#process.load("ACFrameworkModules.MultiCandidateSelector.MultiCandidateSelector_cfi")

### jet energy correction
process.load("ACFrameworkModules.FinalTreeFiller.JetEnergyCorrection_cff")
### MET corrections (produces branch 'pfType1CorrectedMet' which contains type I (L1FastL2L3-L1) corrected pfMET for ak5PFJets with pt>10 GeV)
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
if isData:
    jetCorrectionProducer = process.ak5PFJetsL1L2L3Residual
    process.FinalTreeFiller.pfJets = cms.InputTag('ak5PFJetsL1L2L3Residual')
    process.pfJetMETcorr.jetCorrLabel = cms.string('ak5PFL1FastL2L3Residual')
else:
    jetCorrectionProducer = process.ak5PFJetsL1L2L3
    process.FinalTreeFiller.pfJets = cms.InputTag('ak5PFJetsL1L2L3')

if rerunTauID:
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    tauIDSeq = cms.Sequence(process.PFTau)
    process.FinalTreeFiller.pfTauDiscriminatorPattern = cms.InputTag('hpsPFTauDiscriminationBy','',processName)
else:
    tauIDSeq = cms.Sequence()

#process.debugOutput = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('debugOutput.root'),
#)
#process.out_step = cms.EndPath(process.debugOutput)

#ignore filter
process.ignorePath = cms.Path(
    process.generatorSequence
    *tauIDSeq
    *jetCorrectionProducer
    *process.producePFMETCorrections
    *cms.ignore(process.PrimVtxSelector)
    #*cms.ignore(process.HLTSelector)
    *cms.ignore(process.InputTrackSelector)
    *process.ignoreThreeProngInputSelector
    *cms.ignore(process.KinematicTauProducer)
    *process.matchingSeq
    #*process.MultiCandidateSelector
    *process.FinalTreeFiller
)
#cumulative
process.cumulativePath = cms.Path(
    process.generatorSequence
    *tauIDSeq
    *jetCorrectionProducer
    *process.producePFMETCorrections
    *process.PrimVtxSelector
    #*process.HLTSelector
    *process.InputTrackSelector
    *process.ThreeProngInputSelector
    *process.KinematicTauProducer
    #*process.KinematicTauBasicProducer
    #*process.KinematicTauSkim
    *process.matchingSeq
    #*process.MultiCandidateSelector
    *process.FinalTreeFiller
)

#if no schedule is provided all pathes are executed
if ignoreFilter:
	process.schedule = cms.Schedule(process.ignorePath)
else:
	process.schedule = cms.Schedule(process.cumulativePath)

#process.schedule.append(process.out_step)
