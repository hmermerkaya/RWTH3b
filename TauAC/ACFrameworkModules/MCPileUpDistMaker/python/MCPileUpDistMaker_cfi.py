import FWCore.ParameterSet.Config as cms
import ACFrameworkModules.MCPileUpDistMaker.PileUpDistributions_cfi as MCPileUpDistributions
import os

baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")

MCPileUpDistMaker = cms.EDAnalyzer('MCPileUpDistMaker',
    pileUpDistributions = MCPileUpDistributions.PileUpDistributions,
    outputFile = cms.string("file://"+baseDir+"/data/MCPileUpDistributions.root")
)
