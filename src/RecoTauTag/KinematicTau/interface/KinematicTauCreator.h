// -*- C++ -*-
//
// Package:    KinematicTauCreator
// Class:      KinematicTauCreator
// 
/**
 * This is a pure abstract class providing the interface to the KinematicFit. Use derived classes e.g. ThreeProngTauCreator to implement the create() function.
 *
 * @author Lars Perchalla, Philip Sauerland in 2010
 */
//
// $Id: KinematicTauCreator.h,v 1.15 2010/08/13 12:00:08 perchall Exp $
//
//
//
// Original Author:  Lars Perchalla, Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id: KinematicTauCreator.h,v 1.15 2010/08/13 12:00:08 perchall Exp $
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
	/**
	 default constructor
	 */
    KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder);
	/**
	 constructor with user defined fit parameters
	 */
    KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder, const edm::ParameterSet& cfg);
    virtual ~KinematicTauCreator();

	/**
	 main function of the KinematicTau algorithm. It will call the Fitter and provide the fitted results. It returns true if the fit was successfull.
	 */
    virtual int create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks) = 0;

	/**
	 visible tau constructed from refitted tracks excluding neutrals (the stored primary vertex might be rotated)
	 */
    reco::PFTau getPFTau() const;
	/**
	 the full tau constructed from refitted tracks including neutrals (the stored primary vertex might be rotated)
	 */
	reco::PFTau getKinematicTau() const;
	/**
	 returns vector of all refitted charged particles involved in the fit
	 */
    std::vector<math::XYZTLorentzVector> getRefittedChargedDaughters() const;
	/**
	 returns vector of all refitted neutral particles involved in the fit
	 */
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
	 If the primary vertex was modified by a call of create() this function returns a fake vertex at the modified position with the initial errors. 
	 Otherwise an invalid reco::Vertex is returned.
	 */
	reco::Vertex getModifiedPrimaryVertex() const;
	
protected:
	/**
	 each call of create() will update these members
	 */
    KinematicConstrainedVertexFitter *kcvFitter_;
    RefCountedKinematicTree kinTree_;
	/**
	 ref to original tracks used by the fit
	 */
    std::vector<reco::TrackRef> selectedTracks_;
	TransientTrackBuilder transTrackBuilder_;
	/**
	 Fake vertex. The primary vertex can be rotated within its errors around the secondary vertex. Errors are copied from the original vertex as they have not been modified.
	 */
	reco::Vertex modifiedPV_;
};
