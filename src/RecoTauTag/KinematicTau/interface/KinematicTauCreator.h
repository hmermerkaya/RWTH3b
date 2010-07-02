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
// $Id: KinematicTauCreator.h,v 1.12 2010/06/10 15:43:30 perchall Exp $
//
//
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.12 2010/06/10 15:43:30 perchall Exp $
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
	 visible tau constructed from refitted tracks excluding neutrals (the stored primary vertex might be rotated)
	 */
    reco::PFTau getPFTau() const;
	/**
	 the full tau constructed from refitted tracks including neutrals (the stored primary vertex might be rotated)
	 */
	reco::PFTau getKinematicTau() const;
    std::vector<math::XYZTLorentzVector> getRefittedChargedDaughters() const;
    std::vector<math::XYZTLorentzVector> getRefittedNeutralDaughters() const;

	/**
	 ref to original tracks used by the fit
	 */
    std::vector<reco::TrackRef> getSelectedTracks() const;
	
	/**
	 for debug issues only
	 */
    RefCountedKinematicTree getKinematicTree() const;
	
	/**
	 for debug issues only
	 */
	KinematicConstrainedVertexFitter * getFitter() const {return kcvFitter_;}
	
	/**
	 If the primary vertex was modified by call of create() this function returns a fake vertex at the modified position with the initial errors. 
	 Otherwise an invalid reco::Vertex is returned.
	 */
	reco::Vertex getModifiedPrimaryVertex() const;
	
protected:
	/**
	 each call of create() will update these members
	 */
    KinematicConstrainedVertexFitter *kcvFitter_;
    RefCountedKinematicTree kinTree_;
    std::vector<reco::TrackRef> selectedTracks_;
	TransientTrackBuilder transTrackBuilder_;
	/**
	 fake vertex. the primary vertex can be rotated within its errors around the secondary vertex. errors are copied from original vertex as they have not been modified.
	 */
	reco::Vertex modifiedPV_;
};
