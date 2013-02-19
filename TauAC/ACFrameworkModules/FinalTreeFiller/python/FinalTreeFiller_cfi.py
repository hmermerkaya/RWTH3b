import FWCore.ParameterSet.Config as cms

FinalTreeFiller = cms.EDAnalyzer('FinalTreeFiller',
    genSignal = cms.InputTag('GenSelector','genSignalDecay'),
    genSignalRef = cms.InputTag('GenSelector','genSignalDecayRef'), #same length/content as genSignalDecay, but originals not copied objects
    chargedTauDaughterMatchMap = cms.InputTag('chargedTauDaughterTruth'),
	beamSpot = cms.InputTag('offlineBeamSpot'), #derived from multi-track fit (>>100 tracks)
	primVtx = cms.InputTag('offlinePrimaryVertices'), #this is the PV from the standard reco
#	reducedPrimVtx = cms.InputTag('ThreeProngInputSelectorStep2','primVtx'), #this is the vertex obtained by ThreeProngInputSelector by ignoring the tau tracks (this one is not rotated!!!)
	reducedPrimVtx = cms.InputTag('reducedPrimaryVertices'), #this is the FULL vertex collection obtained by ThreeProngInputSelector by ignoring the tau tracks (this one is not rotated!!! possibly not all of them entered the fit)
    pileupInfo = cms.InputTag('addPileupInfo'),
    triggerResults = cms.InputTag('TriggerResults','','HLT'),
    triggerEvent = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
    #tracks = cms.InputTag('MultiCandidateSelector','generalTracks'),
    tracks = cms.InputTag('generalTracks'),
    minTrackPt = cms.untracked.double(1.),   #ignore tracks below this pt threshold
    muons = cms.InputTag('muons'),
    electrons = cms.InputTag('gsfElectrons'),
    kinematicTaus = cms.InputTag('KinematicTauProducer','SelectedKinematicDecays'),
    pfMET = cms.InputTag('pfMet'),
    pfType1CorrectedMet = cms.InputTag('pfType1CorrectedMet'),
    tcMET = cms.InputTag('tcMet'),
    pfJets = cms.InputTag('ak5PFJets'),
    minJetPt = cms.untracked.double(10.),   #ignore jets below this pt threshold
    pfTaus = cms.InputTag('hpsPFTauProducer'), #the most recently produced product
    minTauPt = cms.untracked.double(10.),   #ignore pftaus below this pt threshold (default is 0., should be defined according to InputTrackSelector)
    pfTauDiscriminatorPattern = cms.InputTag('hpsPFTauDiscriminationBy','','RECO'),# just a pattern (all branches of type PFTauDiscriminator that match this pattern are kept)
    #to be done: replace this string by a scan of all executed modules and test whether it is a filter or not
	flags = cms.vstring("GenSelector","PrimVtxSelector","InputTrackSelector","ThreeProngInputSelectorStep1","ThreeProngInputSelectorStep2","KinematicTauProducer"),
    pileUpDistributionFileMC = cms.untracked.string(""),
    pileUpDistributionHistMC = cms.untracked.string("pileup"),
    pileUpDistributionFileData = cms.untracked.string(""),
    pileUpDistributionHistData = cms.untracked.string("pileup"),
    era = cms.untracked.string("era"), #datataking period that defines the choice of object IDs and so forth
)

chargedTauDaughterTruth = cms.EDProducer("MCTruthDeltaRMatcherNew",#keep both collection as small as possible to prevent fake matches 
    src = cms.InputTag("KinematicTauProducer","selectedTauDaughters"),
    distMin = cms.double(0.01),	#dR=0.1
    matched = FinalTreeFiller.genSignal,
    matchPDGId = cms.vint32(11,13,211),
)

matchingSeq = cms.Sequence(chargedTauDaughterTruth)
