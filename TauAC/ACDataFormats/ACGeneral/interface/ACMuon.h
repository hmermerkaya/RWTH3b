#ifndef ACMuon_h
#define ACMuon_h
/**
 Class to store offline muons according to framework object (reco::Muon).
 More sophisticated constructor in correpsonding conversion class.
 
 Isolation parameters refer to pfIsolationR04 as recommended for 2011, 2012.

 @author Lars Perchalla
 @date 2012
 */

#include <map>

#include "ACCandidate.h"
#include "ACRef.h"

class ACMuon : public ACCandidate {
public:
    virtual std::string classname() const { return "ACMuon"; }
    ACMuon();
    ACMuon(const TLorentzVector& p4, int charge);
    virtual ~ACMuon() {};

    const TVector3 & vertex() const { return vertex_; }
    const ACTrackRef & trackRef() const { return trackRef_; }
    
    bool isGlobalMuonPromptTight() const {return isGlobalMuonPromptTight_;}
    bool isRecommendedMuon() const {return isRecommendedMuon_;}
    
    float sumChargedHadronPt() const {return sumChargedHadronPt_;}
    float sumChargedParticlePt() const {return sumChargedParticlePt_;}
    float sumNeutralHadronEt() const {return sumNeutralHadronEt_;}
    float sumPhotonEt() const {return sumPhotonEt_;}
    float sumPUPt() const {return sumPUPt_;}
//    float emEt() const {return emEt_;}
//    float emVetoEt() const {return emVetoEt_;}
//    float hadEt() const {return hadEt_;}
//    float hadVetoEt() const {return hadVetoEt_;}
//    float hoEt() const {return hoEt_;}
//    float hoVetoEt() const {return hoVetoEt_;}
//    int nJets() const {return nJets_;}
//    int nTracks() const {return nTracks_;}
//    float sumPt() const {return sumPt_;}
//    float trackerVetoPt() const {return trackerVetoPt_;}
    
    /// x coordinate of vertex position
    double vx() const;
    /// y coordinate of vertex position
    double vy() const;
    /// z coordinate of vertex position
    double vz() const;

    /// calculate the isolation according to the POG recommendation for 2011, 2012
    double isolation() const;
    /// calculate the dBeta corrected isolation according to the POG recommendation for 2011, 2012
    double isolationDBeta() const;
    /// whether the muon passes the recommended PF isolation for 2011, 2012
    bool isIsolatedLoose() const;
    /// whether the muon passes the recommended PF isolation for 2011, 2012
    bool isIsolatedTight() const;

    ClassDef(ACMuon,1);

protected:
    /// vertex position
    TVector3 vertex_;
    /// reference to initial track
    ACTrackRef trackRef_;
    /// reference to generator muon (if available)
    ACGenParticleRef genRef_;
    
    /// whether the muon passes muon::GlobalMuonPromptTight (part of the recommendation)
    bool isGlobalMuonPromptTight_;
    /// whether the muon passes the POG recommendation for 2011/2012 (includes isGlobalMuonPromptTight_)
    bool isRecommendedMuon_;
    
    ///muon isolation: sum-pt of charged Hadron
    float sumChargedHadronPt_;
    ///muon isolation: sum-pt of charged Particles(inludes e/mu)
    float sumChargedParticlePt_;
    ///muon isolation: sum pt of neutral hadrons
    float sumNeutralHadronEt_;
    ///muon isolation: sum pt of PF photons
    float sumPhotonEt_;
    ///muon isolation: sum pt of charged Particles not from PV (for Pu corrections)
    float sumPUPt_;

//deprecated:
//    ///muon isolation: ECAL sum-Et
//    float emEt_;
// 	///muon isolation: ECAL sum-et in the veto region in r-phi
//    float emVetoEt_;
// 	///muon isolation: hcal sum-Et
//    float hadEt_;
// 	///muon isolation: hcal sum-et in the veto region in r-phi
//    float hadVetoEt_;
// 	///muon isolation: ho sum-Et
//    float hoEt_;
// 	///muon isolation: ho sum-et in the veto region in r-phi
//    float hoVetoEt_;
// 	///muon isolation: number of jets in the cone
//    int nJets_;
// 	///muon isolation: number of tracks in the cone (excluding veto region)
//    int nTracks_;
// 	///muon isolation: sum-pt of tracks
//    float sumPt_;
// 	///muon isolation: (sum-)pt inside the veto region in r-phi
//    float trackerVetoPt_;
};

typedef ACMuon* pACMuon;

#endif
