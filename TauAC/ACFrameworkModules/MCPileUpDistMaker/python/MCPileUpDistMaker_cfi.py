import FWCore.ParameterSet.Config as cms
import ACFrameworkModules.MCPileUpDistMaker.PileUpDistributions_cfi as MCPileUpDistributions

MCPileUpDistMaker = cms.EDAnalyzer('MCPileUpDistMaker',
    pileUpDistributions = MCPileUpDistributions.PileUpDistributions,
    outputFile = cms.string("file://../FinalTreeFiller/MCPileUpDistributions.root")
)
