import os
import hashlib
import subprocess

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

baseDir = os.path.relpath(os.environ.get('CMSSW_BASE')+"/src")

def createDataPileUpFile(jsonFile, pileUpReweightingType):
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()
    thisHash = hashlib.md5(lumiString)
    #print 'Require pileup file with hash ', thisHash.hexdigest()
    dataPileUpFilename = '/data/DataPileUp_'+thisHash.hexdigest()+'_'+pileUpReweightingType+'.root'
    if not os.path.exists(baseDir+dataPileUpFilename):
        print 'Creating new "', pileUpReweightingType, '" data pileup file for json file "',jsonFile,'"'
        subprocess.call(['pileupCalc.py', '-i', jsonFile, '--inputLumiJSON', '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt', '--calcMode', pileUpReweightingType, '--maxPileupBin', '50', baseDir+dataPileUpFilename], shell=False)
    return dataPileUpFilename
