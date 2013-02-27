#ifndef ParticleBuilder_h
#define ParticleBuilder_h

#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "TString.h"
#include "TVector3.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class ParticleBuilder {
 public:
  enum CMSSWPerigee{aCurv=0,aTheta,aPhi,aTip,aLip};

  ParticleBuilder(){};
  ~ParticleBuilder(){};

  static TrackParticle CreateTrackParticle(const reco::TrackRef &track,  edm::ESHandle<TransientTrackBuilder>  &transTrackBuilder, const GlobalPoint p,bool fromPerigee=true, bool useTrackHelixPropogation=true);
  static reco::Vertex  GetVertex(LorentzVectorParticle p);

 private:
  static TMatrixT<double> ConvertCMSSWTrackParToSFTrackPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ConvertCMSSWTrackPerigeeToSFTrackPar(TMatrixT<double> &inpar);

};
#endif


