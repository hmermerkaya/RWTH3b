//
// Original Author:  Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.7 2010/01/20 17:41:06 perchall Exp $
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

	/**
	 visible tau constructed from refitted tracks
	 */
    reco::PFTau getPFTau();
    std::vector<math::XYZTLorentzVector> getRefittedChargedHadrons();

	/**
	 ref to original tracks used by the fit
	 */
    std::vector<reco::TrackRef> getSelectedTracks();
	
	/**
	 for debug issues only
	 */
    RefCountedKinematicTree getKinematicTree();
	
	/**
	 for debug issues only
	 */
	KinematicConstrainedVertexFitter * getFitter(){return kcvFitter_;}
	
protected:
    KinematicConstrainedVertexFitter *kcvFitter_;
    RefCountedKinematicTree kinTree_;
    std::vector<reco::TrackRef> selectedTracks_;
	TransientTrackBuilder transTrackBuilder_;
};
