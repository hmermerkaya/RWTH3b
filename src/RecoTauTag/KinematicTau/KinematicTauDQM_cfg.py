

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
process.GlobalTag.globaltag = 'START44_V9B::All'

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
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6EB56796-403D-E111-B113-001EC9D83141.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6EB56796-403D-E111-B113-001EC9D83141.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/FE9A1BA3-413D-E111-AA25-90E6BA19A22E.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/ECD8D599-413D-E111-9428-485B39800BD7.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/ACA5AF97-413D-E111-8C68-20CF305B0519.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/76360A97-413D-E111-B614-20CF3027A5ED.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/A6ABB3AC-413D-E111-99CF-001EC9D446AD.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/1A99E739-423D-E111-AC54-00261834B5AF.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/3C29673A-423D-E111-A2A4-E0CB4E29C514.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/2244D93A-423D-E111-9967-E0CB4E553642.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/2A76E3FF-423D-E111-BAC4-20CF305616E1.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6A8F05FD-423D-E111-BA9D-E0CB4E553642.root')
                            )


numberOfEvents = -1
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

process.endjob_step = cms.Path(process.endOfProcess)
process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.postProcessorKinematicFitValidation*process.dqmSaver)
process.schedule = cms.Schedule(process.KinFitSkim,process.endjob_step)
