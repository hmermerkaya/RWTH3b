#setup every parameter that depends on a certain user

def LoadUserParameters():
    username = 'perchalla'
    outputPath = '/user/'+username+'/output/TauAC/CMSSW_4_4_0/debug/'
    jsonFile = "./Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt"

    inputPath = 'Summer11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S4_START42_V11-v1/'
    jobName = 'tauSMBrFromZ'

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

    return [inputPath, outputPath, jobName, isData, jsonFile]
