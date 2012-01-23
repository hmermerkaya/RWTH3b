import os
import hashlib

#parse lumi regions
from FWCore.PythonUtilities.LumiList import LumiList

def createDataPileUpFile(jsonFile):
    jsonList = LumiList(jsonFile)
    lumiString = jsonList.getCMSSWString()
    thisHash = hashlib.md5(lumiString)
    dataPileUpFilename = 'DataPileUp_'+thisHash.hexdigest()+'.root'
    if not os.path.exists(dataPileUpFilename):
        print 'Creating new data pileup file for json file "',jsonFile,'"'
        os.system('pileupCalc.py -i '+jsonFile+' --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/pileup_JSON_2011_4_2_validation.txt --calcMode observed '+dataPileUpFilename)
    return dataPileUpFilename
