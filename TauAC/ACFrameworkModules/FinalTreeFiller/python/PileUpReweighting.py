import os
import hashlib
import subprocess

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

def createDataPileUpFile(jsonFile):
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()
    thisHash = hashlib.md5(lumiString)
    dataPileUpFilename = 'DataPileUp_'+thisHash.hexdigest()
    dataPileUpFilenameTrue = dataPileUpFilename+'_true.root'
    dataPileUpFilenameObs = dataPileUpFilename+'_obs.root'
    if not os.path.exists(dataPileUpFilenameObs):
        print 'Creating new "observed" data pileup file for json file "',jsonFile,'"'
        subprocess.call(['pileupCalc.py', '-i', jsonFile, '--inputLumiJSON', '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt', '--calcMode', 'observed', '--maxPileupBin', '50', dataPileUpFilenameObs], shell=False)
    if not os.path.exists(dataPileUpFilenameTrue):
        print 'Creating new "true" data pileup file for json file "',jsonFile,'"'
        subprocess.call(['pileupCalc.py', '-i', jsonFile, '--inputLumiJSON', '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt', '--calcMode', 'true', '--maxPileupBin', '50', dataPileUpFilenameTrue], shell=False)
    return [dataPileUpFilenameTrue, dataPileUpFilenameObs] 
