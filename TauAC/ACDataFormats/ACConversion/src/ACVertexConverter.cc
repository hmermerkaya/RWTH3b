#include "../interface/ACVertexConverter.h"

ACVertexConverter::ACVertexConverter(const reco::Vertex& vertex) {
    position_ = TVector3(vertex.x(), vertex.y(), vertex.z());
    valid_ = vertex.isValid();
    error_.ResizeTo(TMatrixDSym(3));
    for (int i=0; i!=3; i++)    for (int j=0; j!=3; j++) error_(i,j) = vertex.covariance(i,j);
    chi2_ = vertex.chi2();
    ndof_ = vertex.ndof();
    trackSize_ = vertex.tracksSize();
    struct VertexHigherPtSquared vtxPtCalculation;
    ptSum_ = vtxPtCalculation.sumPtSquared(vertex);
}
