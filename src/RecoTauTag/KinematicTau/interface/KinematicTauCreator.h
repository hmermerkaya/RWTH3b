//
// Original Author:  Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id$
//
//


#include "FWCore/ParameterSet/interface/ParameterSet.h"



class KinematicTauCreator
{
public:
    explicit KinematicTauCreator(const edm::ParameterSet& cfg=0);
    ~KinematicTauCreator();
        
protected:
    virtual int create(const reco::Vertex& primvtx_, const std::vector<reco::TrackRef>& inputTracks);
    reco::PFTau getPFTau();
    std::vector<math::XYZTLorentzVector> getRefittedChargedHadrons();
    RefCountedKinematicTree getKinematicTree();
    int iterations();
    float csum();

    KinematicConstrainedVertexFitter kcvFitter;
    RefCountedKinematicTree kinTree;
    int iterations;
    float csum;
};
