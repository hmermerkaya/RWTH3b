import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

process = cms.Process("AOD")

process.load('Configuration/StandardSequences/Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# global tag
#process.GlobalTag.globaltag = 'START44_V9B::All'
process.GlobalTag.globaltag = 'START52_V9::All'
#process.GlobalTag.globaltag = 'FT_R_53_V6C::All'



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
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/56D5F502-C295-E111-A388-003048D46004.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/88EA3847-8E96-E111-894F-00E0817917C7.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/78A70594-7796-E111-9B82-003048D4772E.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/2883A54B-3E97-E111-BA1A-001A647894E8.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/C09EAADE-3D97-E111-B32E-001E67396DF1.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/D2E4D132-3D97-E111-88B2-003048673FC0.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/32696749-3E97-E111-BB99-001A647894A0.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/02F49D49-3E97-E111-9AF6-003048670A0C.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/521DD948-9F96-E111-B544-003048D45FD8.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/EA57A1DE-CD96-E111-8D2F-0025B3E063A8.root",
    "file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/760F45FA-3D97-E111-A86A-002590200974.root",
    )
                            )
numberOfEvents = 1000

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(numberOfEvents)
        )

process.options = cms.untracked.PSet(
        Rethrow = cms.untracked.vstring('ProductNotFound')
        )


process.load("RecoTauTag.KinematicTau.KinematicFitSequences_cff")
process.load("RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi")
#process.KinematicTauProducer.useTrackHelixFit = cms.untracked.bool(True)
                                      

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
process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.schedule = cms.Schedule(process.KinFitSkim,process.endjob_step)

#sequences for Kinematic Fit 
#KinematicFitSequence         = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer)
#KinematicFitSequencewithSkim = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer*KinematicTauSkim)
#KinematicFitSequencewithDQM  = cms.Sequence(ThreeProngInputSelectorStep1*KinematicTauProducer*KinematicTauAnalyzer)
