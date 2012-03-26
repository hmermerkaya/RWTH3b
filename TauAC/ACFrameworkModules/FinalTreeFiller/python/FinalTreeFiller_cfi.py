import FWCore.ParameterSet.Config as cms

FinalTreeFiller = cms.EDAnalyzer('FinalTreeFiller',
    genSignal = cms.InputTag('GenSelector','genSignalDecay'),
    genSignalRef = cms.InputTag('GenSelector','genSignalDecayRef'), #same length/content as genSignalDecay, but originals not copied objects
    chargedTauDaughterMatchMap = cms.InputTag('chargedTauDaughterTruth'),
	primVtx = cms.InputTag('offlinePrimaryVertices'), #this is the PV from the standard reco
#	reducedPrimVtx = cms.InputTag('ThreeProngInputSelectorStep2','primVtx'), #this is the vertex obtained by ThreeProngInputSelector by ignoring the tau tracks (this one is not rotated!!!)
	reducedPrimVtx = cms.InputTag('reducedPrimaryVertices'), #this is the FULL vertex collection obtained by ThreeProngInputSelector by ignoring the tau tracks (this one is not rotated!!! possibly not all of them entered the fit)
    pileupInfo = cms.InputTag('addPileupInfo'),
    triggerResults = cms.InputTag('TriggerResults','','HLT'),
    tracks = cms.InputTag('MultiCandidateSelector','generalTracks'),
    muons = cms.InputTag('muons'),
    electrons = cms.InputTag('gsfElectrons'),
    kinematicTaus = cms.InputTag('KinematicTauProducer','SelectedKinematicDecays'),
    pfMET = cms.InputTag('pfMet'),
    tcMET = cms.InputTag('tcMet'),
    pfJets = cms.InputTag('MultiCandidateSelector','ak5PFJets'),
    pfTaus = cms.InputTag('hpsPFTauProducer'),
    pfTauDiscriminators = cms.VInputTag(
        cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
        cms.InputTag('hpsPFTauDiscriminationByLooseElectronRejection'),
        cms.InputTag('hpsPFTauDiscriminationByLooseIsolation'),
        cms.InputTag('hpsPFTauDiscriminationByLooseMuonRejection'),
        cms.InputTag('hpsPFTauDiscriminationByMediumElectronRejection'),
        cms.InputTag('hpsPFTauDiscriminationByMediumIsolation'),
        cms.InputTag('hpsPFTauDiscriminationByTightElectronRejection'),
        cms.InputTag('hpsPFTauDiscriminationByTightIsolation'),
        cms.InputTag('hpsPFTauDiscriminationByTightMuonRejection'),
        cms.InputTag('hpsPFTauDiscriminationByVLooseIsolation'),
    ),
    #to be done: replace this string by a scan of all executed modules and test whether it is a filter or not
	flags = cms.vstring("GenSelector","PrimVtxSelector","InputTrackSelector","ThreeProngInputSelectorStep1","ThreeProngInputSelectorStep2","KinematicTauProducer"),
    pileUpDistributionFileMC = cms.untracked.string(""),
    pileUpDistributionHistMC = cms.untracked.string("pileup"),
    pileUpDistributionFileData = cms.untracked.string(""),
    pileUpDistributionHistData = cms.untracked.string("pileup"),
)

chargedTauDaughterTruth = cms.EDProducer("MCTruthDeltaRMatcherNew",#keep both collection as small as possible to prevent fake matches 
    src = cms.InputTag("KinematicTauProducer","selectedTauDaughters"),
    distMin = cms.double(0.01),	#dR=0.1
    matched = FinalTreeFiller.genSignal,
    matchPDGId = cms.vint32(11,13,211),
)

matchingSeq = cms.Sequence(chargedTauDaughterTruth)
