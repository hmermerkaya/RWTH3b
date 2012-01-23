import FWCore.ParameterSet.Config as cms

process = cms.Process("MCPileUpDistMaker")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
                        input = cms.untracked.int32(0)
                    )

process.load("ACFrameworkModules.MCPileUpDistMaker.MCPileUpDistMaker_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = process.MCPileUpDistMaker.outputFile,
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.MCPileUpDistMaker
)

process.schedule = cms.Schedule(process.p)


