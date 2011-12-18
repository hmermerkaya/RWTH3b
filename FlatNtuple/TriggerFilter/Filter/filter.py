import FWCore.ParameterSet.Config as cms

process = cms.Process("Filter")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.GlobalTag.globaltag = 'START3X_V26::All'
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")


numberOfEvents = 10

inPath1 = '/user/cherepanov/data_CMSSW_3_8_7/9d2fd75bcd33af5022152f60e3b68073/'  # Z tau tau number oif files = 403
f1=open('/user/cherepanov/data_CMSSW_3_8_7/filelist.dat')

listOfFiles=[]


nFiles1 = 491
for nf in range(1,nFiles1):
    string=f1.readline()
    listOfFiles.append('file://'+inPath1+string[:-1] )
#    print nf

#print listOfFiles
      
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#    'dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Spring10/W1Jets_Pt0to100-alpgen/GEN-SIM-RECO/START3X_V26_S09-v1/0025/423A407A-9048-DF11-9B59-002618943937.root'),
                            'dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms//store/mc/Fall10/Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/GEN-SIM-RECO/START38_V12-v1/0018/BC50C6FA-A9E2-DF11-8A75-0018FE284C0E.root'),
#        listOfFiles),
	noEventSort = cms.untracked.bool(True), #Events are processed in order of run number and luminosity block number always. If this parameter is true, then within a lumino	duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)
    

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numberOfEvents)
)

process.load("TriggerFltr.TrigFilter.TrigFilter_cfi")



process.p = cms.Path(process.TrigFilter)
