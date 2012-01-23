import FWCore.ParameterSet.Config as cms

from FWCore.MessageLogger.MessageLogger_cfi import *

#print 'destinations:',MessageLogger.destinations
#print 'statistics:',MessageLogger.statistics
#print 'categories:',MessageLogger.categories
#print 'debugModules:',MessageLogger.debugModules
#print 'cerr:',MessageLogger.cerr
#print 'cout:',MessageLogger.cout
#print 'cerr_stats',MessageLogger.cerr_stats

useDEBUG = False
enableStatistics = False

severity = 'INFO'
if useDEBUG:
  severity = 'DEBUG'


#MessageLogger.destinations = cms.untracked.vstring(
#'warnings','errors','infos','debugs',#standard destinations printed to file
#'cout','cerr'#standard destinations printed to screen. own extensions create a new file
#)
if enableStatistics:
  MessageLogger.statistics.append('test_statistics')#default is cerr_stats

#see infos on debug modules at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
#to enable debugging messages use:
#scram b USER_CXXFLAGS="-DEDM_ML_DEBUG"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MessageLogger.debugModules = cms.untracked.vstring('HLTSelector','KinematicTauProducer','KinematicTauAdvancedProducer','InputTrackSelector','ThreeProngInputSelector','GenSelector','GenJetSelector')#name of the modules which may deliver messages
#categories for all severity levels are defined in the log message: LogDebug("category")<<... One needs to append a category to listen to it
MessageLogger.categories.append('HLTSelector')
MessageLogger.categories.append('KinematicTauCreator')
MessageLogger.categories.append('ThreeProngInputSelector')
MessageLogger.categories.append('InputTrackSelector')
MessageLogger.categories.append('KinematicTauAdvancedProducer')
MessageLogger.categories.append('ThreeProngTauCreator')
MessageLogger.categories.append('KinematicConstrainedVertexFitter')
MessageLogger.categories.append('KinematicConstrainedVertexUpdator')
MessageLogger.categories.append('ConstrainedTreeBuilder')
MessageLogger.categories.append('MCCollectionNotFound')#some modules requiring MC information throw this kind of messages when run on data

MessageLogger.cerr = cms.untracked.PSet(
	noTimeStamps = cms.untracked.bool(True),
  threshold = cms.untracked.string(severity),#es gibt 4 severity-stufen in folgender reihenfolge: "DEBUG", "INFO"(default), "WARNING", "ERROR"
  FwkReport = cms.untracked.PSet(#manipuliert ausgabe von main_input:source, welches evt vorliegt
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(0)
  )
)

if useDEBUG:
    print '-->To enable debugging compile with: scram b USER_CXXFLAGS="-DEDM_ML_DEBUG"'
    MessageLogger.cerr = cms.untracked.PSet(#use cerr/cout for printing to screen, debug for creating a new file
		  threshold = cms.untracked.string(severity),
		  FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
		  FwkSummary = cms.untracked.PSet(limit = cms.untracked.int32(0)),
      #INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),#disable INFO in general
      #DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),#disable DEBUG in general
      #HLTSelector = cms.untracked.PSet(limit = cms.untracked.int32(0)),#listen only to specific module(s)
	)
#print '-->FinalTreeFiller.MessageLogger_cfi::Warning: Some messages will be partially suppressed if like "ThreeProngTauCreator::kinematicRefit: ERROR! Tree is not valid."'
#MessageLogger.cerr.ThreeProngTauCreator = cms.untracked.PSet(limit = cms.untracked.int32(1))#reduce reoccurence (after limit is reached reports are shown exponentially)
MessageLogger.cerr.MCCollectionNotFound = cms.untracked.PSet(limit = cms.untracked.int32(1))

if enableStatistics:
  MessageLogger.test_statistics = cms.untracked.PSet(
   output = cms.untracked.string('cerr'),
   threshold = cms.untracked.string(severity),
   FwkJob = cms.untracked.PSet(limit = cms.untracked.int32(0)),
   FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),                                                     
   FwkSummary = cms.untracked.PSet(limit = cms.untracked.int32(0)),
   #WARNING = cms.untracked.PSet(limit = cms.untracked.int32(0)),
   #INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
   #only listen to specific debug module
   #DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
   #ThreeProngInputSelector = cms.untracked.PSet(limit = cms.untracked.int32(-1))
  )
