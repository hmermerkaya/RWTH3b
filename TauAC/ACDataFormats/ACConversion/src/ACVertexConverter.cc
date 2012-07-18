#include "../interface/ACVertexConverter.h"

ACVertexConverter::ACVertexConverter(const reco::Vertex& vertex, const TrackConversionLog * conversionLogTrack) {
    position_ = TVector3(vertex.x(), vertex.y(), vertex.z());
    valid_ = vertex.isValid();
    error_.ResizeTo(TMatrixDSym(3));
    for (int i=0; i!=3; i++)    for (int j=0; j!=3; j++) error_(i,j) = vertex.covariance(i,j);
    chi2_ = vertex.chi2();
    ndof_ = vertex.ndof();
    trackSize_ = vertex.tracksSize();
    struct VertexHigherPtSquared vtxPtCalculation;
    ptSum_ = vtxPtCalculation.sumPtSquared(vertex);
    
    if (!conversionLogTrack) return; //might be missing due to different track collection
    for (std::vector<reco::TrackBaseRef>::const_iterator trk = vertex.tracks_begin(); trk != vertex.tracks_end(); ++trk) {
        const ACTrackRef * trkRef = 0;
        if (trk->isAvailable()) {
            try {
                // This function throws an exception if the type passed does not match the concrete reference type
                reco::TrackRef concrete = trk->castTo<reco::TrackRef>();
                if (concrete.isAvailable()) {
                    trkRef = conversionLogTrack->getConverted(concrete);
                }
            } catch (edm::errors::ErrorCodes error) {
                std::cout<<"ACVertexConverter::ACVertexConverter: Catched an error "<<error<<" (expecting "<<edm::errors::InvalidReference<<")"<<std::endl;
            }
        }
        // reference might be missing (e.g. due to candidate thresholds)
        if (trkRef) {
            tracks_.push_back(*trkRef);
            //std::cout<<"Found valid ref to track with pt "<<(*trk)->pt()<<std::endl;
        }
    }
    //std::cout<<"ACVertexConverter::ACVertexConverter: Found "<<tracks_.size()<<"/"<<trackSize()<<" valid refs to tracks."<<std::endl;
}
