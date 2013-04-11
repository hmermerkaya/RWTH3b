
import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

process = cms.Process("AOD")

process.load('Configuration/StandardSequences/Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# global tag
#process.GlobalTag.globaltag = 'START52_V9::All'
process.GlobalTag.globaltag = 'GR_R_52_V10::All'
#process.GlobalTag.globaltag = 'GR_R_52_V4::All'


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('HLTrigger.Configuration.HLT_GRun_cff')

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

######################################################

process.MessageLogger.categories.append('InputTrackSelector')
process.MessageLogger.categories.append('KinematicTau')
process.MessageLogger.categories.append('KinematicTauSkim')
process.MessageLogger.categories.append('ThreeProngInputSelector_Step2')
process.MessageLogger.debugModules = cms.untracked.vstring('KinematicTau', 'KinematicTauSkim','InputTrackSelector','ThreeProngInputSelector_Step2')
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    InputTrackSelector = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    KinematicTau = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    KinematicTauSkim = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    ThreeProngInputSelector_Step2 = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/home/home2/institut_3b/nugent/2883A54B-3E97-E111-BA1A-001A647894E8.root',
    'file:/home/home2/institut_3b/nugent/28A87D76-5EC2-E111-923C-001D09F23174.root')
                            
                            
                            )

numberOfEvents = 60
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
    )

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
    )


process.load("RecoTauTag.KinematicTau.KinematicFitSequences_cff")
process.load("RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi")


process.load('Configuration.StandardSequences.EDMtoMEAtJobEnd_cff')
process.load("Validation.Configuration.postValidation_cff")

process.schedule = cms.Schedule()

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtJobEnd_cff')
process.load('RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi')
process.load('RecoTauTag.KinematicTau.Tau_JAKID_Filter_cfi')
process.load('DQMServices.Components.MEtoEDMConverter_cfi')

process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.dqmSaver.workflow = "/KinematicFitSequencewithDQM/VAL/RECO"
#process.DQMStore.verbose=1


process.KinematicTauAnalyzer.doFakeRate  = cms.bool(True)
process.KinematicTauAnalyzer.doDQM = cms.bool(False)
                                      

process.endjob_step = cms.Path(process.endOfProcess)
#process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.KinFitSkim  = cms.Path(process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.schedule = cms.Schedule(process.KinFitSkim,process.endjob_step)
