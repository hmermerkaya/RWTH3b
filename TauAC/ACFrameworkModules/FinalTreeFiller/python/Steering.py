import os
import sys
import FWCore.ParameterSet.VarParsing as VarParsing

#setup every parameter that depends on a certain user

baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")

def registerCommandlineArguments():
    #see default arguments at
    #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#Passing_Command_Line_Arguments_T
    
    #print "arguments are:\n", sys.argv
    options = VarParsing.VarParsing ('analysis')
    options.register ('dataset',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "unique name of the dataset to run on (impacts global tag, generator selection, ...)")
    options.register ('era',
                  2011, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,         # string, int, or float
                  "datataking period (impacts e.g. the object IDs according to particular recommendation)")
    options.register ('runFromCrab',
                      0, # default value
                      VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                      VarParsing.VarParsing.varType.int,         # string, int, or float
                      "When running via crab use runFromCrab=1 (default is 0).")
    options.register ('lumiJSON',
                      "", # default value
                      VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                      VarParsing.VarParsing.varType.string,         # string, int, or float
                      "JSON to filter lumis (data, if run locally)")

    # get and parse the command line arguments
    options.parseArguments()
    
    #validate options
    if options.era != 2011 and options.era != 2012:
        print "ERROR! Unexpected commandline option era =", options.era
        sys.exit(404)
    
    print "--> CMSSW will run with the following command line arguments:"
    print " "*3, "dataset =", options.dataset
    print " "*3, "era =", options.era
    print " "*3, "lumiJSON =", options.lumiJSON if options.lumiJSON else "DEFAULT"
    
    return options

def LoadUserParameters(options):
    username = 'perchalla'
    outputPath = '/user/'+username+'/output/analysis/CMSSW_4_4_2/development/'
    if options.era == 2012:
        outputPath = '/user/'+username+'/output/TauAC/CMSSW_5_3_8/development/'
    
    jobName = options.dataset
    if not jobName:
        print "ERROR! No datset selected! Pass the dataset definition as command line argument by dataset=<myDataset>"
        sys.exit(404)

    jsonFile = options.lumiJSON #only needed for local lumi masking (will be ignored in crab mode)
    pileUpDistributionMC = ''
    if options.era==2011:
        pileUpDistributionMC = 'MCPileUpDistMaker/Fall2011'
        if not jsonFile:
            jsonFile = baseDir+"/ACFrameworkModules/FinalTreeFiller/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt"
    elif options.era==2012:
        pileUpDistributionMC = ''
        if not jsonFile:
            jsonFile = baseDir+"/ACFrameworkModules/FinalTreeFiller/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt"
    else:
        print "ERROR! Unknown datataking period", options.era
        sys.exit(404)

    inputPath = ''
    globaltag = ''
    isData = False #data or MC

    if '2011' in jobName or '2012' in jobName:
        isData = True
    if options.era == 2011:
        if isData: globaltag = 'GR_R_44_V15::All'
        else:      globaltag = 'START44_V12::All'
    elif options.era == 2012:
        if isData:
            #2012 A and B: data re-reco with CMSSW_5_3_2_patch4
            if '2012A' in jobName or '2012B' in jobName: globaltag = 'FT_53_V6_AN3::All'
            #2012 C v1: data re-reco with CMSSW_5_3_2_patch4
            if '2012C1' in jobName: globaltag = 'FT_53_V10_AN3::All'
            #Global tags for prompt reco data analysis (starting with 2012)
            if '2012C2' in jobName: globaltag = 'GR_P_V41_AN3::All'
            if '2012D' in jobName: globaltag = 'GR_P_V42_AN3::All'
        else: globaltag = 'START53_V7G::All' #2012 MC production: Summer12_DR53X (START53_V7A used during production)
    if not globaltag:
        print "ERROR! No GlobalTag specified!"
        sys.exit(404)
    else:
        print "--> Analysing", "data" if isData else "MC", "using the GlobalTag", globaltag

    #the following is only required for local execution
    if jobName == 'DYTauTau':
        ##/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM, globaltag START42_V11::All (from evtDump, dbs query fails here), HLT: /online/collisions/2011/5e32/v6.2/HLT/V4
        #inputPath = 'Summer11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S4_START42_V11-v1/'
        #globaltag = 'START44_V7::All'
        ##/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall11-PU_S6-START44_V5-v1/AODSIM, produced with GlobalTag START44_V5::All
        inputPath = 'Fall11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S6-START44_V5-v1/'
    elif jobName == 'QCD30to80':
        ##/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/Summer11-PU_S4_START42_V11-v1/AODSIM, globaltag START42_V11::All (also from evtDump, dbs query fails here too), HLT: /online/collisions/2011/5e32/v6.2/HLT/V4
        #inputPath = 'Summer11/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/AODSIM'
        #globaltag = 'START44_V7::All'
        ##/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/Fall11-PU_S6-START44_V5-v1/AODSIM
        inputPath = 'Fall11/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia/AODSIM/PU_S6-START44_V5-v1/'
    elif jobName == 'QCD80to170':
        ##/QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'QCD170to250':
        ##/QCD_Pt-170to250_EMEnriched_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'QCD250to350':
        ##/QCD_Pt-250to350_EMEnriched_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'QCD350toInf':
        ##/QCD_Pt-350_EMEnriched_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'WJetsToLNu':
        ##/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6-START44_V5-v1/AODSIM
        ##or /WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'DYJetsToLL':
        ##/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'WW':
        ##/WW_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'WZ':
        ##/WZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'ZZ':
        ##/ZZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'TTJets':
        ##/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'TauNov2011A':
        ##/Tau/Run2011A-08Nov2011-v1/AOD
        inputPath = 'Run2011A/Tau/AOD/08Nov2011-v1/'
        isData = True
    elif jobName == 'TauNov2011B':
        ##/Tau/Run2011B-19Nov2011-v1/AOD
        inputPath = ''
        isData = True
    elif jobName == 'TauPlusXNov2011A':
        ##/TauPlusX/Run2011A-08Nov2011-v1/AOD
        inputPath = 'Run2011A/TauPlusX/AOD/08Nov2011-v1/'
        isData = True
    elif jobName == 'TauPlusXNov2011B':
        ##/TauPlusX/Run2011A-19Nov2011-v1/AOD
        inputPath = ''
        isData = True
    elif jobName == 'MinimumBias':
        ##/MinimumBias/Run2011A-PromptReco-v6/AOD   globaltag='GR_P_V22::All'
        inputPath = 'Run2011A/MinimumBias/AOD/PromptReco-v6/'
        isData = True
    elif jobName == 'GFHTauTauM115':
        ##/GluGluToHToTauTau_M-115_7TeV-powheg-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'GFHTauTauM125':
        ##/GluGluToHToTauTau_M-125_7TeV-powheg-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = ''
    elif jobName == 'GFHTauTauM135':
        ##/GluGluToHToTauTau_M-135_7TeV-powheg-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = 'Fall11/GluGluToHToTauTau_M-135_7TeV-powheg-pythia6/AODSIM/PU_S6_START42_V14B-v1/'
    elif jobName == 'SingleMuNov2011A':
        ##/SingleMu/Run2011A-08Nov2011-v1/AOD
        inputPath = 'Run2011A/SingleMu/AOD/08Nov2011-v1/'
        isData = True
    elif jobName == 'SingleMuNov2011B':
        ##/SingleMu/Run2011B-19Nov2011-v1/AOD
        inputPath = ''
        isData = True
    elif jobName == 'Parking2012B':
        inputPath = 'Run2012B/ParkingMonitor/AOD/'
        isData = True
    elif jobName == 'TauPlusX2012A':
        inputPath = 'Run2012A/TauPlusX/AOD/'
        isData = True
    elif jobName == 'TTJetsSemiL':
        ##/GluGluToHToTauTau_M-135_7TeV-powheg-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM
        inputPath = 'Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/'


    disk = 'disk1'
    #disk = 'user'
    
    mainputPath = '/'+disk+'/'+username+'/store/'
    if not isData:
        mainputPath += 'mc/'
    else:
        mainputPath += 'data/'
    inputPath = mainputPath + inputPath

    return [inputPath, outputPath, jobName, globaltag, isData, jsonFile, pileUpDistributionMC]


def returnOptionsString():
    """ format the options to be passed on the command line"""
    outputString = ""
    for optionsName, optionValue in options.__dict__['_singletons'].iteritems():
        outputString += " %s=%s" % (optionsName, optionValue)
    
    for optionsName, optionValues in options.__dict__['_lists'].iteritems():
        for anOption in optionValues:
            outputString += " %s=%s" % (optionsName, anOption)
    return outputString
