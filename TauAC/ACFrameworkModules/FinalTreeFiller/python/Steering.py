import os
#setup every parameter that depends on a certain user

baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")

def LoadUserParameters():
    username = 'perchalla'
    outputPath = '/user/'+username+'/output/analysis/CMSSW_4_4_2/PVBugFix/'

#    jsonFile = baseDir+"/ACFrameworkModules/FinalTreeFiller/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt"
    jsonFile = baseDir+"/ACFrameworkModules/FinalTreeFiller/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt"

    jobName = 'TauNov2011'

    inputPath = ''
    globaltag = ''
#    pileUpDistributionMC = 'MCPileUpDistMaker/Summer11_PU_S4'
    pileUpDistributionMC = 'MCPileUpDistMaker/Fall2011'

    if jobName == 'tauSMBrFromZ':
        ##/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM, globaltag START42_V11::All (from evtDump, dbs query fails here), HLT: /online/collisions/2011/5e32/v6.2/HLT/V4
        #inputPath = 'Summer11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S4_START42_V11-v1/'
        #globaltag = 'START44_V7::All'
        ##/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall11-PU_S6-START44_V5-v1/AODSIM, produced with GlobalTag START44_V5::All
        inputPath = 'Fall11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S6-START44_V5-v1/'
        globaltag = 'START44_V12::All'
    elif jobName == 'QCD':
        ##/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/Summer11-PU_S4_START42_V11-v1/AODSIM, globaltag START42_V11::All (also from evtDump, dbs query fails here too), HLT: /online/collisions/2011/5e32/v6.2/HLT/V4
        #inputPath = 'Summer11/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/AODSIM'
        #globaltag = 'START44_V7::All'
        ##/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/Fall11-PU_S6-START44_V5-v1/AODSIM
        inputPath = 'Fall11/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/AODSIM/PU_S6-START44_V5-v1/'
        globaltag = 'START44_V12::All'
    elif jobName == 'TauNov2011':
        ##/Tau/Run2011A-08Nov2011-v1/AOD
        inputPath = 'Run2011A/Tau/AOD/08Nov2011-v1/'
        globaltag = 'GR_R_44_V14::All'
    elif jobName == 'TauPlusXNov2011':
        ##/TauPlusX/Run2011A-08Nov2011-v1/AOD
        inputPath = 'Run2011A/TauPlusX/AOD/08Nov2011-v1/'
        globaltag = 'GR_R_44_V14::All'
    elif jobName == 'MinimumBias':
        ##/MinimumBias/Run2011A-PromptReco-v6/AOD   globaltag='GR_P_V22::All'
        inputPath = 'Run2011A/MinimumBias/AOD/PromptReco-v6/'
        globaltag = 'GR_P_V27::All'
    
    isData = False #data or MC
    if inputPath.find('Run20', 0, 5) != -1:
        isData = True

    disk = 'disk1'
    #disk = 'user'
    
    mainputPath = '/'+disk+'/'+username+'/store/'
    if not isData:
        mainputPath += 'mc/'
    else:
        mainputPath += 'data/'
    inputPath = mainputPath + inputPath

    return [inputPath, outputPath, jobName, globaltag, isData, jsonFile, pileUpDistributionMC]
