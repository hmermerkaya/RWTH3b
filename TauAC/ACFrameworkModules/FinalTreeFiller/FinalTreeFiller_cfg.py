import FWCore.ParameterSet.Config as cms
import os
import sys
import ACFrameworkModules.Common.Tools as Tools
import ACFrameworkModules.FinalTreeFiller.Steering as Steering
import ACFrameworkModules.FinalTreeFiller.PileUpReweighting as PileUpReweighting

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

#relative path to $CMSSW_BASE/src. needed for crab and to run from arbitrary paths
baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")
release = os.environ.get('CMSSW_VERSION')

processName = "FinalTreeFiller"
process = cms.Process(processName)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("ACFrameworkModules.Common.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

options = Steering.registerCommandlineArguments()
###############
#steering parameters#
inputPath, outputPath, jobName, globaltag, isData, jsonFile, pileUpDistributionMC = Steering.LoadUserParameters(options)
numberOfEvents = -1
ignoreFilter = False # but always ignore TauMotherProducer until it is ready
printEvents = 0	#print generator event
prepareEvtDisplay = False # whether or not to store KinematicFit's output to edm branches for easy drawing in cmsShow
triggerTag = "HLT"
doPileUpReweighting = True
pileUpReweightingType = "true" #(true or observed) choose according to https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities
minTau = 1 #minimum of selected taus
minJetPt = 10. #ignore pfjets below this pt threshold (in GeV)
minTauPt = 10. #ignore pftaus below this pt threshold (in GeV)
minTrackPt = 1. #ignore tracks below this pt threshold (in GeV)
rerunTauID = True #whether or not to rerun the PFTau sequence
if options.era != 2011: rerunTauID = False #in CMSSW_5_3_x stick to the release for now
###############
print ' '*3, 'Run', release, 'with the following cuts: minTau =', minTau, ', minJetPt =', minJetPt, ', minTauPt =', minTauPt, ', minTrackPt =', minTrackPt
print ' '*3, 'The TauID is', 'reexecuted' if rerunTauID else 'taken from the release.'

if options.era == 2011:
    if release.lstrip('CMSSW_')[0] != '4':
        print 'ERROR! Release', release, 'is not supposed to run on era', options.era
        sys.exit(1)
elif options.era == 2012:
    if release.lstrip('CMSSW_')[0] != '5':
        print 'ERROR! Release', release, 'is not supposed to run on era', options.era
        sys.exit(1)
else:
    print 'ERROR! Unexpected era', options.era
    sys.exit(1)


#summaries paths and modules (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMPathsAndTriggerBits#Summary_report)
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) #not needed as EventLogger stores the path summary

if not options.runFromCrab:
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

### select runs from JSON file
if isData and not options.runFromCrab:
    print 'Discriminate lumisections from', jsonFile
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()#convert into compact format needed by CMSSW
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(lumiString.split(','))

process.GlobalTag.globaltag = globaltag
#print 'Using the GlobalTag ',process.GlobalTag.globaltag

motherPdgID = 0
if isData:
    decayType = 'data'
else:
    decayType, motherPdgID = Tools.parseDecayType(jobName)

numberOfEventsStr = str(numberOfEvents)
if numberOfEvents == -1: numberOfEventsStr = 'all'

print '--> Analysing',numberOfEventsStr,'events in <'+inputPath+'>.'
print ' '*3, 'Detected',len(sourceString) ,'files of decayType = '+decayType+'(',motherPdgID,').'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('file://'+outputPath+jobName+'_'+numberOfEventsStr+'.root'), #+date.today().strftime('%Y%m%d')+'.root'),
    closeFileFast = cms.untracked.bool(True)
)
if options.runFromCrab:
    process.TFileService.fileName = cms.string(jobName+'_'+numberOfEventsStr+'.root')

### cleaning according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCollisionsDataAnalysis#Recipes_to_get_started
process.noscraping = cms.EDFilter("FilterOutScraping",
	applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24), 
    maxd0 = cms.double(2) 
)
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.cleaning = cms.Sequence(process.noscraping*process.primaryVertexFilter*process.HBHENoiseFilter)


###skim workflow
process.load("ACFrameworkModules.FinalTreeFiller.TriggerResultsFilter_cfi")
process.triggerSelection.hltResults = cms.InputTag("TriggerResults","", triggerTag)
if options.era == 2011:
    process.triggerSelection.triggerConditions = cms.vstring('HLT_IsoMu15_*LooseIsoPFTau20_v*')
    #or in 2011 data allow also for IsoMu12_LooseIsoPFTau10, IsoMu15_LooseIsoPFTau15 and IsoMu15_eta2p1_LooseIsoPFTau20?
if not isData and options.era == 2012:
    process.triggerSelection.triggerConditions = cms.vstring('HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*')
if not ignoreFilter: print " "*3, "Apply HLT filter on path(s):", ", ".join(process.triggerSelection.triggerConditions)

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
process.load("RecoTauTag.KinematicTau.PrepareEvtDisplay_cfi")

### event storage
process.load("ACFrameworkModules.FinalTreeFiller.FinalTreeFiller_cfi")
process.FinalTreeFiller.triggerResults = cms.InputTag("TriggerResults","",triggerTag)
process.FinalTreeFiller.decayType = cms.untracked.string(decayType)
process.FinalTreeFiller.minTauPt = cms.untracked.double(minTauPt)
process.FinalTreeFiller.minJetPt = cms.untracked.double(minJetPt)
process.FinalTreeFiller.minTrackPt = cms.untracked.double(minTrackPt)
process.FinalTreeFiller.era = cms.untracked.string(str(options.era))
process.load("ACFrameworkModules.EventLogger.EventLogger_cfi")

if doPileUpReweighting and not isData:
    #dataPileUpFilename = baseDir + PileUpReweighting.createDataPileUpFile(jsonFile, pileUpReweightingType)
    #process.FinalTreeFiller.pileUpDistributionFileMC = cms.untracked.string(baseDir+"/data/MCPileUpDistributions.root")
    #process.FinalTreeFiller.pileUpDistributionHistMC = cms.untracked.string(pileUpDistributionMC)
    #process.FinalTreeFiller.pileUpDistributionFileData = cms.untracked.string(dataPileUpFilename)
    #process.FinalTreeFiller.pileUpDistributionHistData = cms.untracked.string("pileup")
    process.FinalTreeFiller.pileUpDistributionHistMC = cms.untracked.string("pileup")
    process.FinalTreeFiller.pileUpDistributionHistData = cms.untracked.string("pileup")
    process.EventLogger.pileUpDistributionHistMC = cms.untracked.string("pileup")
    process.EventLogger.pileUpDistributionHistData = cms.untracked.string("pileup")
    pileupMCHist = ""
    pileupDataHist = ""
    if options.era == 2011:
        pileupMCHist = "/data/MC_Fall11_PU_S6-500bins.root"
        pileupDataHist = "/data/Data_Pileup_2011_HCP-500bins.root"
    if options.era == 2012:
        pileupMCHist = "/data/MC_Summer12_PU_S10-600bins.root"
        pileupDataHist = "/data/Data_Pileup_2012_Moriond-600bins.root"
    process.FinalTreeFiller.pileUpDistributionFileMC = cms.untracked.string(baseDir+pileupMCHist)
    process.FinalTreeFiller.pileUpDistributionFileData = cms.untracked.string(baseDir+pileupDataHist)
    process.EventLogger.pileUpDistributionFileMC = cms.untracked.string(baseDir+pileupMCHist)
    process.EventLogger.pileUpDistributionFileData = cms.untracked.string(baseDir+pileupDataHist)
    print " "*3, "Apply pileup reweighting from histograms:", pileupMCHist, "(MC),", pileupDataHist, "(data)"


### candidate selectors (deprecated as violating refs)
#process.load("ACFrameworkModules.MultiCandidateSelector.MultiCandidateSelector_cfi")

### jet energy correction
process.load("ACFrameworkModules.FinalTreeFiller.JetEnergyCorrection_cff")
### Type-1 MET corrections (for Jet Energy Scale) (produces branch 'pfType1CorrectedMet' which contains type I (L1FastL2L3-L1) corrected pfMET for ak5PFJets with pt>10 GeV)
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
if isData:
    jetCorrectionProducer = process.ak5PFJetsL1L2L3Residual
    process.FinalTreeFiller.pfJets = cms.InputTag('ak5PFJetsL1L2L3Residual')
    process.pfJetMETcorr.jetCorrLabel = cms.string('ak5PFL1FastL2L3Residual')
else:
    jetCorrectionProducer = process.ak5PFJetsL1L2L3
    process.FinalTreeFiller.pfJets = cms.InputTag('ak5PFJetsL1L2L3')
#Type-II MET Correction (for Unclustered Energy) is NOT RECOMMENDED for PFMet --> skip it
#Type-0 MET Correction with PFCandidate (for Pile-up) (hack to work on RECO instead of PAT). This will add the type-0 correction to 'pfType1CorrectedMet'
if options.era == '2012':
    process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
    process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
    process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    )
    process.METCorrections = cms.Sequence(process.type0PFMEtCorrection*process.producePFMETCorrections)
else:
    process.METCorrections = cms.Sequence(process.producePFMETCorrections)
#MET x/y Shift Correction (for phi modulation) - NOT YET APPLIED



if rerunTauID:
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    tauIDSeq = cms.Sequence(process.PFTau)
    process.FinalTreeFiller.pfTauDiscriminatorPattern = cms.InputTag('hpsPFTauDiscriminationBy','',processName)
else:
    tauIDSeq = cms.Sequence()

if prepareEvtDisplay:
    process.debugOutput = cms.OutputModule("PoolOutputModule",
        outputCommands = cms.untracked.vstring('keep *'),
        fileName = cms.untracked.string('debugOutput.root'),
    )
    process.out_step = cms.EndPath(process.debugOutput)

#ignore filter
process.ignorePath = cms.Path(
    process.generatorSequence
    *tauIDSeq
    *jetCorrectionProducer
    *process.METCorrections
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
    process.cleaning
    *process.triggerSelection
    *process.generatorSequence
    *tauIDSeq
    *jetCorrectionProducer
    *process.METCorrections
    *process.PrimVtxSelector
    *process.InputTrackSelector
    *process.ThreeProngInputSelector
    *process.KinematicTauProducer
    #*process.KinematicTauBasicProducer
    #*process.KinematicTauSkim
    *process.matchingSeq
    #*process.MultiCandidateSelector
    *process.FinalTreeFiller
)

scanner = Tools.Scanner()
#if no schedule is provided all pathes are executed
if ignoreFilter:
    process.schedule = cms.Schedule(process.ignorePath)
    process.EventLogger.pathName = process.ignorePath.label() # monitor only the main path
    process.cumulativePath.visit(scanner)
else:
    process.schedule = cms.Schedule(process.cumulativePath)
    process.EventLogger.pathName = process.cumulativePath.label() # monitor only the main path
    process.cumulativePath.visit(scanner)

#scanner.printModules()
#print scanner.filterLabels()
#get rid of all those EDAnalyzers and EDProducers (especially those from rerunning the TauID)
process.EventLogger.modules = scanner.filterLabels() # only monitor EDFilter (all other types of modules may not veto events)

if prepareEvtDisplay:
    process.ignorePath *= process.PrepareEvtDisplay
    process.cumulativePath *= process.PrepareEvtDisplay
    process.schedule.append(process.out_step)

process.log_step = cms.EndPath(process.EventLogger)
process.schedule.append(process.log_step)
