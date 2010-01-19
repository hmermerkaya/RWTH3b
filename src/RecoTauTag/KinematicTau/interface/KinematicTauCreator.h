//
// Original Author:  Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.4 2010/01/15 17:15:46 perchall Exp $
//
//


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


class KinematicTauCreator
{
public:
    KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder);
    KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder, const edm::ParameterSet& cfg);
    virtual ~KinematicTauCreator();

    virtual int create(const reco::Vertex& primvtx_, const std::vector<reco::TrackRef>& inputTracks) = 0;
    reco::PFTau getPFTau();
    std::vector<math::XYZTLorentzVector> getRefittedChargedHadrons();
    RefCountedKinematicTree getKinematicTree();
    int iterations();
    float csum();
	
protected:
    KinematicConstrainedVertexFitter *kcvFitter_;
    RefCountedKinematicTree kinTree_;
    int iterations_;
    float csum_;
	TransientTrackBuilder transTrackBuilder_;
};
