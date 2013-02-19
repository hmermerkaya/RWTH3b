import os
import fnmatch
import FWCore.ParameterSet.Config as cms

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
        if 'FromZ' in jobName or 'DYTauTau' in jobName:
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


class Scanner(object):
    """Class to scan a sequence/path and store its modules. (See visit() function in /FWCore/ParameterSet/python/SequenceTypes.py for enter and leave).
        Usage: scanner = Scanner()
               mysequence.visit(scanner)
               scanner.printModules()
        """
    def __init__(self):
        self.modules = []
    def enter(self,visitee):
        """to be called from the visit function"""
        self.modules.append(visitee)
    def modules(self):
        """return the stored objects"""
        return self.modules
    def filterLabels(self):
        """return the labels of stored EDFilter objects"""
        labels = [module.label() for module in self.modules if type(module) is cms.EDFilter]
        return labels
    def filters(self):
        """return the stored objects of kind EDFilter"""
        modules = [module for module in self.modules if type(module) is cms.EDFilter]
        return modules
    def analyzers(self):
        """return the stored objects of kind EDAnalyzer"""
        modules = [module for module in self.modules if type(module) is cms.EDAnalyzer]
        return modules
    def producers(self):
        """return the stored objects of kind EDProducer"""
        modules = [module for module in self.modules if type(module) is cms.EDProducer]
        return modules
    def printModules(self):
        """print the stored objects"""
        print "Scanner stores", len(self.filters()),"EDFilter", len(self.producers()),"EDProducer", len(self.analyzers()),"EDAnalyzer"
        for module in self.modules:
            print "module ", module, "is filter?", "YES" if type(module) is cms.EDFilter else "NO"
            #if type(module) in [cms.EDFilter, cms.EDProducer, cms.EDAnalyzer]:
            #    print "module ", module.label(), "is of type", type(module)
    def leave(self, visitee):
        """to be called from the visit function"""
        pass
