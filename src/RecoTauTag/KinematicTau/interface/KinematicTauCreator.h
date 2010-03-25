// -*- C++ -*-
//
// Package:    KinematicTauCreator
// Class:      KinematicTauCreator
// 
/**
 
 Description: pure abstract class. use e.g. ThreeProngTauCreator to implement create()
 
 Implementation:
 <Notes on implementation>
 */
//
// $Id: KinematicTauCreator.h,v 1.10 2010/01/27 14:28:53 perchall Exp $
//
//
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.10 2010/01/27 14:28:53 perchall Exp $
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

    virtual int create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks) = 0;

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
	
	/**
	 If the primary vertex was modified by call of create() this function returns a fake vertex at the modified position. 
	 Otherwise an invalid reco::Vertex is returned.
	 */
	reco::Vertex getModifiedPrimaryVertex();
	
protected:
	/**
	 each call of create() will update these members
	 */
    KinematicConstrainedVertexFitter *kcvFitter_;
    RefCountedKinematicTree kinTree_;
    std::vector<reco::TrackRef> selectedTracks_;
	TransientTrackBuilder transTrackBuilder_;
	/**
	 fake vertex. the primary vertex can be rotated within its errors around the secondary vertex.
	 */
	reco::Vertex modifiedPV_;
};
