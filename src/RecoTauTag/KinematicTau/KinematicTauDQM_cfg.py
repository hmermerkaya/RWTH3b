

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
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6EB56796-403D-E111-B113-001EC9D83141.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6EB56796-403D-E111-B113-001EC9D83141.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/FE9A1BA3-413D-E111-AA25-90E6BA19A22E.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/ECD8D599-413D-E111-9428-485B39800BD7.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/ACA5AF97-413D-E111-8C68-20CF305B0519.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/76360A97-413D-E111-B614-20CF3027A5ED.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/A6ABB3AC-413D-E111-99CF-001EC9D446AD.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/1A99E739-423D-E111-AC54-00261834B5AF.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/3C29673A-423D-E111-A2A4-E0CB4E29C514.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/2244D93A-423D-E111-9967-E0CB4E553642.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/2A76E3FF-423D-E111-BAC4-20CF305616E1.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/6A8F05FD-423D-E111-BA9D-E0CB4E553642.root')
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/0C7A62FF-842D-E111-94F8-002618943916.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/82920D11-A82D-E111-96ED-002618943896.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/522DE5D9-A62D-E111-9C28-001A92810A94.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/BC2B6F4C-B02D-E111-B784-0026189438CE.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/241AC25E-9A2D-E111-BD64-0018F3D09620.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/BE9BA0C5-932D-E111-B119-0026189438CF.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/C2CBF935-912D-E111-BAD4-003048678F74.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/408BB648-AC2D-E111-B1B8-0026189438AD.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/5819235F-922D-E111-B361-0018F3D096DE.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/485AFFBF-A12D-E111-B433-0018F3D0960E.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/8C2B8DD7-AA2D-E111-9AB3-003048FFCB96.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/74354302-A62D-E111-ACBE-0018F3D096DA.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/907D4682-962D-E111-B1C7-002618FDA287.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/FCEF8FB3-AB2D-E111-8305-002618FDA237.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Fall11/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0001/54A0D4EF-B32D-E111-8555-001A92811724.root')
                           
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/56D5F502-C295-E111-A388-003048D46004.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/88EA3847-8E96-E111-894F-00E0817917C7.root',
    #'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/9626DDBC-FC96-E111-9902-984BE1089EB8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/A88FB007-6996-E111-AD99-003048D477C0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/143B1001-EA95-E111-B59C-001E67397D55.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/CAFE2D50-7B96-E111-A080-003048D4602E.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/5A8C276A-C796-E111-9D70-0025B3E06698.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/4A92A241-6196-E111-AC2D-001E67397003.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/D84355E3-A496-E111-B856-003048673E7A.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/447CB04E-9496-E111-BAF8-00E0817918A7.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/1E329673-4E97-E111-A3BA-002590200A6C.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/CA2F29F4-A096-E111-B982-003048D47796.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/1890E748-3E97-E111-9C53-001E67396CFC.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/9C89B183-3E97-E111-9B63-003048670BCC.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/C6EC5782-8A96-E111-A4CF-003048673EFE.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/5C3DF315-CF96-E111-9323-0025B3E05BF4.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/FEF8EAF4-E096-E111-8A40-003048D477A0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/C0B17BAD-CC96-E111-BF74-003048D4772C.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/2CA03D61-7A96-E111-9D72-003048D462AE.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/F4FFE93C-DE96-E111-8825-003048D477BE.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/5EBE22A1-BE96-E111-A7CA-003048D45FCA.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/B4D432C2-4E96-E111-8E94-001E6739665D.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/EEC15947-9496-E111-93C9-001E67398C5A.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/4CAAF0F5-7696-E111-893C-003048D46006.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/2E676562-4296-E111-AE30-0025B3E05D0A.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/BCC2931C-9896-E111-BBAC-0025902008B8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/FCA63B0D-C296-E111-B5D9-00E08179173F.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/AAEEFCF5-A496-E111-8ECA-00E081791779.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/A425BF2E-F196-E111-A86E-003048D47A0E.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/2655F9F9-7B96-E111-937B-003048D47766.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/8CDCEBD4-C596-E111-89A3-001E67396568.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/EAB97421-ED95-E111-B6D9-001E67397CC9.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/78A70594-7796-E111-9B82-003048D4772E.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/2883A54B-3E97-E111-BA1A-001A647894E8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/C09EAADE-3D97-E111-B32E-001E67396DF1.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/D2E4D132-3D97-E111-88B2-003048673FC0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/32696749-3E97-E111-BB99-001A647894A0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/02F49D49-3E97-E111-9AF6-003048670A0C.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/521DD948-9F96-E111-B544-003048D45FD8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/EA57A1DE-CD96-E111-8D2F-0025B3E063A8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/760F45FA-3D97-E111-A86A-002590200974.root')

                            
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
#process.KinFitSkim  = cms.Path(process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
##process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.schedule = cms.Schedule(process.KinFitSkim,process.endjob_step)
