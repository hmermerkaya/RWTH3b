#include "../interface/ACMuonConverter.h"

ACMuonConverter::ACMuonConverter(const reco::Muon& particle, bool isRecommendedMuon, const TrackConversionLog * conversionLogTrack, ACGenParticleRef genRef) {
    name_ = "muon";
    charge_ = particle.charge();
    p4_ = TLorentzVector(particle.p4().Px(), particle.p4().Py(), particle.p4().Pz(), particle.p4().E());
    vertex_ = TVector3(particle.vx(), particle.vy(), particle.vz());
    //store vertexCovariance!!!
    genRef_ = genRef;
    /// to be ignored:    ACGenDecayRef genDecayRef_
    
    isGlobalMuonPromptTight_ = muon::isGoodMuon(particle, muon::GlobalMuonPromptTight);
    isRecommendedMuon_ = isRecommendedMuon;

    if (particle.isPFIsolationValid()) {
        //const reco::MuonIsolation  iso03 = particle.pfIsolationR03();
        const reco::MuonPFIsolation  iso04 = particle.pfIsolationR04(); // recommended for 2011
        
        sumChargedHadronPt_ = iso04.sumChargedHadronPt;
        sumChargedParticlePt_ = iso04.sumChargedParticlePt;
        sumNeutralHadronEt_ = iso04.sumNeutralHadronEt;
        sumPhotonEt_ = iso04.sumPhotonEt;
        sumPUPt_ = iso04.sumPUPt;
    } else {
        sumChargedHadronPt_ = -1.;
        sumChargedParticlePt_ = -1.;
        sumNeutralHadronEt_ = -1.;
        sumPhotonEt_ = -1.;
        sumPUPt_ = -1.;
    }

//deprecated
//    if (particle.isIsolationValid()) {
//        const reco::MuonIsolation iso03 = particle.isolationR03();
//
//        emEt_ = iso03.emEt();
//        emVetoEt_ = iso03.emVetoEt();
//        hadEt_ = iso03.hadEt();
//        hadVetoEt_ = iso03.hadVetoEt();
//        hoEt_ = iso03.hoEt();
//        hoVetoEt_ = iso03.hoVetoEt();
//        nJets_ = iso03.nJets();
//        nTracks_ = iso03.nTracks();
//        sumPt_ = iso03.sumPt();
//        trackerVetoPt_ = iso03.trackerVetoPt();
//    } else {
//        emEt_ = -1.;
//        emVetoEt_ = -1.;
//        hadEt_ = -1.;
//        hadVetoEt_ = -1.;
//        hoEt_ = -1.;
//        hoVetoEt_ = -1.;
//        nJets_ = -1;
//        nTracks_ = -1;
//        sumPt_ = -1.;
//        trackerVetoPt_ = -1.;
//    }
    
    trackRef_ = ACTrackRef();
    if (!conversionLogTrack) return; //might be missing due to different track collection
    const ACTrackRef * trkRef = 0;
    if (particle.track().isAvailable()) {
        trkRef = conversionLogTrack->getConverted(particle.track());
    }
    // reference might be missing (e.g. due to candidate thresholds)
    if (trkRef) {
        trackRef_ = *trkRef;
    }

}
