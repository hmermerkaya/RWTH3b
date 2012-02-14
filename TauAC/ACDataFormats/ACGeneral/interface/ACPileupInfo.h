#ifndef ACPileupInfo_h
#define ACPileupInfo_h
/**
 Container to store pileup information. One-to-one copy of the corresponding CMSSW class /SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h.

 Note: In the Summer11 samples, the "true" number of pileup interactions that were overlaid on a given event is what is reported by getPU_NumInteractions().
 Previously, the number of interactions per event could have been slightly under-reported because the calculation was based on the number of events that left some particle in the tracking volume.
 With the expansion of the acceptance for forward physics simulation in the Summer11 samples, there are many more events with no visible tracks that can be included.
 This means that the length of the vectors computed from the tracking quantities is less than or equal to the number of interactions in the crossing.
 Use the size() function to figure out how many interactions were visible to the tracker.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCommon.h"

class ACPileupInfo : public ACCommon {
public:
    virtual std::string classname() const { return "ACPileupInfo"; }
    ACPileupInfo();
    ACPileupInfo(const int num_PU_vertices, std::vector<float>& zpositions, std::vector<float>& sumpT_lowpT, std::vector<float>& sumpT_highpT, std::vector<int>& ntrks_lowpT, std::vector<int>& ntrks_highpT);
    ACPileupInfo(const int num_PU_vertices, std::vector<float>& zpositions, std::vector<float>& sumpT_lowpT, std::vector<float>& sumpT_highpT, std::vector<int>& ntrks_lowpT, std::vector<int>& ntrks_highpT, int bunchCrossing);
    ACPileupInfo(const int num_PU_vertices, std::vector<float>& zpositions, std::vector<float>& sumpT_lowpT, std::vector<float>& sumpT_highpT, std::vector<int>& ntrks_lowpT, std::vector<int>& ntrks_highpT, int bunchCrossing, float TrueNumInteractions);

//    ACPileupInfo(const int num_PU_vertices,
//                        std::vector<float>& instLumi,
//                        std::vector<edm::EventID>& eventInfo);

    virtual ~ACPileupInfo() {};

    const int getPU_NumInteractions() const { return num_PU_vertices_; }
    const std::vector<float>& getPU_zpositions() const { return zpositions_; }
    const std::vector<float>& getPU_sumpT_lowpT() const { return sumpT_lowpT_; }
    const std::vector<float>& getPU_sumpT_highpT() const { return sumpT_highpT_; }
    const std::vector<int>& getPU_ntrks_lowpT() const { return ntrks_lowpT_; }
    const std::vector<int>& getPU_ntrks_highpT() const { return ntrks_highpT_; }
    const std::vector<float>& getPU_instLumi() const { return instLumi_; }
//    const std::vector<edm::EventID>& getPU_EventID() const { return eventInfo_; }
    const int getBunchCrossing() const {
        if (isnan(bunchCrossing_)) {
            return -1000.0;
        }
        return bunchCrossing_;
    }
    const float getTrueNumInteractions() const { 
        if (isnan(TrueNumInteractions_)) {
            return -1.0;
        }
        return TrueNumInteractions_;
    }
    
    ClassDef(ACPileupInfo,1);

protected:
    // for "standard" pileup: we have MC Truth information for these

    /// the number of pileup interactions that have been added to the event
    int num_PU_vertices_;
    /// the true primary vertex position along the z axis for each added interaction
    std::vector<float> zpositions_;
    /// the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
    std::vector<float> sumpT_lowpT_;
    /// the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
    std::vector<float> sumpT_highpT_;
    /// the number of tracks originating from each interaction, where track pT > low_cut
    std::vector<int> ntrks_lowpT_;
    /// the number of tracks originating from each interaction, where track pT > high_cut
    std::vector<int> ntrks_highpT_;
    /// to which bunch crossing does this interaction belong?0=in-time crossing, negative ones are early, positive ones are late
    int bunchCrossing_;
    ///  true number of interactions for each event
    float TrueNumInteractions_;


    // for DataMixer pileup, we only have raw information:

    std::vector<float> instLumi_;
//    std::vector<edm::EventID> eventInfo_;

};

typedef ACPileupInfo* pACPileupInfo;

#endif
