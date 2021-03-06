#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "TLorentzVector.h"
#include <vector>

namespace {
  struct dictionary {
    SelectedKinematicParticle p1;
    std::vector< SelectedKinematicParticle > vp1;
    std::vector< std::vector< SelectedKinematicParticle > > vvp1;
    edm::Wrapper< std::vector< SelectedKinematicParticle > > wvp1;
    edm::Wrapper< std::vector< std::vector< SelectedKinematicParticle > > > wvvp1;
    
    std::vector<reco::TrackRefVector> vvr2;
    edm::Wrapper<std::vector<reco::TrackRefVector> > wvvr2;
    edm::Wrapper<reco::PFTauRef> wvr3;

    std::vector<reco::Track> vt1;
    edm::Wrapper<std::vector<reco::Track> > wvt1;
    edm::Wrapper<reco::Vertex> wvt2;
    edm::Wrapper<reco::VertexRef> wvt3;    

    std::vector<std::vector< reco::TrackRef > > vvr1;        
    edm::Wrapper<std::vector<std::vector< reco::TrackRef > > > wvvr1;
    std::vector<std::vector<std::vector< reco::TrackRef > > > vvvr1;        
    edm::Wrapper<std::vector<std::vector<std::vector< reco::TrackRef > > > > wvvvr1;
    
    SelectedKinematicDecay d1;
    std::vector< SelectedKinematicDecay > vd1;
    edm::Wrapper< std::vector< SelectedKinematicDecay > > wvd1;
    edm::Wrapper< std::vector<std::vector< SelectedKinematicDecay > > > wvvd1;		
 
    std::map<std::string, bool> map1;
    edm::Wrapper<std::map<std::string, bool> > wmap1;
    std::vector<std::map<std::string, bool> > vmap1;
    edm::Wrapper<std::vector<std::map<std::string, bool> > > wvmap1;

    std::vector<TLorentzVector> vlv;
    edm::Wrapper<std::vector<TLorentzVector> > wvlv;


    std::vector<TVector3> vtdv;
    edm::Wrapper<std::vector<TVector3> > wvtdv;


  };
}

