import FWCore.ParameterSet.Config as cms

#from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *

#details at https://twiki.cern.ch/twiki/bin/view/CMS/TriggerResultsFilter

triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),   # HLT results   - set to empty to ignore HLT
    l1tResults = cms.InputTag( "" ),                            # L1 GT results (gtDigis) - set to empty to ignore L1
    l1tIgnoreMask = cms.bool( False ),                          # use L1 mask
    l1techIgnorePrescales = cms.bool( False ),                  # read L1 technical bits from PSB#9, bypassing the prescales
    daqPartitions = cms.uint32( 1 ),                            # used by the definition of the L1 mask
    throw = cms.bool( False )                                   # True: throw exception on unknown trigger names, False: only report a warning for each pattern in triggerConditions and continue
)

#trigger condition: 2012 data for mu+Tau selection
triggerSelection.triggerConditions = cms.vstring(
    # trigger condition combined by logical OR
    'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*',                  # MuTau selection
    'HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*',
    #'HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*',        # TauTau selection
    #'HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*',
    #'HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*',
)
