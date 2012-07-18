import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *

# this should load the correction service and the corresponding producers
# one producer is e.g. ak5PFJetsL1L2L3 (recommendated for MC) and ak5PFJetsL1L2L3Residual (recommendated for data)
# run one of these to create a new branch which contains the corrected jets
