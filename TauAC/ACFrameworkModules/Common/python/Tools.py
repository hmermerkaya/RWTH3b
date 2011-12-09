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
    if jobName.find('tau3pi')!=-1:
        if jobName.find('tau3piFromVBFH')!=-1:
            decayType = 'VBFH3pr'
        if jobName.find('tau3piFromGFH')!=-1:
            decayType = 'GFH3pr'
        if jobName.find('tau3piFromZ')!=-1:
            decayType = 'Z3pr'
        if jobName.find('tau3piFromW')!=-1:
            decayType = 'W3pr'
    else:  
        if jobName.find('FromZ')!=-1:
            decayType = 'TauPair'
            motherPdgID = 23
        if jobName.find('FromW')!=-1:
            print 'Decay mode not yet covered!'
        #decayType = 'TauPair'
        #motherPdgID = 23
        if jobName.find('FromGFH')!=-1 or jobName.find('FromVBFH')!=-1:
            decayType = 'TauPair'
            motherPdgID = 25

    return [decayType, motherPdgID]
