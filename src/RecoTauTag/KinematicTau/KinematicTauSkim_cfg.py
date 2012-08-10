import FWCore.ParameterSet.Config as cms

process = cms.Process("AODQualityTauSkim")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START36_V9::All'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('InputTrackSelector')
process.MessageLogger.categories.append('KinematicTauProducer')
process.MessageLogger.categories.append('KinematicTauSkim')
process.MessageLogger.debugModules = cms.untracked.vstring('KinematicTauProducer', 'KinematicTauSkim')
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	InputTrackSelector = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
	KinematicTauProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
	KinematicTauSkim = cms.untracked.PSet(limit = cms.untracked.int32(-1))
)

numberOfEvents = -1

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
#	splitLevel = cms.untracked.int32(0),
	outputCommands = process.AODSIMEventContent.outputCommands,
	fileName = cms.untracked.string('QCD_EMEnriched_Pt30to80_AODQualityTauSkim.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('skim_QualityTau')
    )
)
process.output.outputCommands.extend(
	cms.untracked.vstring(
			'keep *_fixedConeHighEffPFTauProducer*_*_*',
			'keep *_fixedConeHighEffPFTauDiscrimination*_*_*'
	)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
		#'dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/perchall/QCD_EMEnriched_Pt30to80/TauSkim/0a829d6435628d6d858898be444d1f11/QCD_EMEnriched_Pt30to80_AODTauSkim_1_1_IUg.root'
		'file:///user/perchalla/data/CMSSW_3_6_2/QCD_EMEnriched_Pt30to80/TauSkim/0a829d6435628d6d858898be444d1f11/QCD_EMEnriched_Pt30to80_AODTauSkim_1_1_IUg.root'
	)
#	,duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

process.load("CommonTools.PrimVtxSelector.PrimVtxSelector_cfi")
process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")
process.load("RecoTauTag.KinematicTau.ThreeProngInputSelector_cff")
process.load("RecoTauTag.KinematicTau.kinematictau_cfi")
process.load("RecoTauTag.KinematicTau.KinematicTauSkim_cfi")

process.skim_QualityTau = cms.Path(process.PrimVtxSelector*process.InputTrackSelector*process.ThreeProngInputSelector*process.KinematicTauProducer*process.KinematicTauSkim)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.skim_QualityTau)
process.schedule.append(process.out_step)
