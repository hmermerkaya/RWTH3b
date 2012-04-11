import os
import fnmatch

def scanDir(dir, pattern='*.root'): #pattern = '*.py' # Can include any UNIX shell-style wildcards
  #scan <dir> recursivly and return any files matching <pattern>
  matches = []
  for root, dirs, files in os.walk(dir):
    for filename in fnmatch.filter(files, pattern):
      matches.append(os.path.join(root, filename))
  
  return matches

def buildPoolSource(files):
  #expects a list of files and inserts the prefix "file://" to prepare a PoolSource
  sourceString = []
  for file in files:
    sourceString.append("file://"+file)
  
  return sourceString

def parseDecayType(jobName):
    decayType = 'unknown'
    motherPdgID = 0
    if 'tau3pi' in jobName:
        if 'tau3piFromVBFH' in jobName:
            decayType = 'VBFH3pr'
        if 'tau3piFromGFH' in jobName:
            decayType = 'GFH3pr'
        if 'tau3piFromZ' in jobName:
            decayType = 'Z3pr'
        if 'tau3piFromW' in jobName:
            decayType = 'W3pr'
    else:  
        if 'FromZ' in jobName:
            decayType = 'TauPair'
            motherPdgID = 23
        if 'FromGFH' in jobName or 'FromVBFH' in jobName or 'GFHTauTau' in jobName:
            decayType = 'TauPair'
            motherPdgID = 25
        if 'FromW' in jobName:
            print 'Decay mode not yet covered!'
        #decayType = 'TauPair'
        #motherPdgID = 23

    #print 'jobname, decayType, motherPdgID: ',jobName, decayType, motherPdgID
    return [decayType, motherPdgID]
