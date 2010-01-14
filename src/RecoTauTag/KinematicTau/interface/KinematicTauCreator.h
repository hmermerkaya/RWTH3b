//
// Original Author:  Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.1 2010/01/14 10:37:17 sauerlan Exp $
//
//


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"


class KinematicTauCreator
{
public:
    KinematicTauCreator();
    KinematicTauCreator(const edm::ParameterSet& cfg);
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
