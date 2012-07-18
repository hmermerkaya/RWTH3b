#include "../interface/ACBeamSpotConverter.h"

ACBeamSpotConverter::ACBeamSpotConverter(const reco::BeamSpot& beamSpot) {
    position_ = TVector3(beamSpot.position().X(), beamSpot.position().Y(), beamSpot.position().Z());

    error_.ResizeTo(TMatrixDSym(7));
    for (int i=0; i!=7; i++) for (int j=0; j!=7; j++) error_(i,j) = beamSpot.covariance(i,j);
    
    sigmaZ_ = beamSpot.sigmaZ();
    beamWidthX_ = beamSpot.BeamWidthX();
    beamWidthY_ = beamSpot.BeamWidthY();
    dxdz_ = beamSpot.dxdz();
    dydz_ = beamSpot.dydz();
    emittanceX_ = beamSpot.emittanceX();
    emittanceY_ = beamSpot.emittanceY();
    betaStar_ = beamSpot.betaStar();
    type_ = beamSpot.type();
}
