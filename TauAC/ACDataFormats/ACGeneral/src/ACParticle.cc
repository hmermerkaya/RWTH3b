#include "../interface/ACParticle.h"

ClassImp(ACParticle)

ACParticle::ACParticle():
ACGenParticle() {
    matrix_.ResizeTo(TMatrixDSym(7));
    genRef_ = ACGenParticleRef();
}
ACParticle::ACParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex, const TMatrixDSym & matrix, ACGenParticleRef genRef, ACGenDecayRef genDecayRef):
ACGenParticle(status, name, pdgId, charge, p4, vertex, genDecayRef) {
    matrix_.ResizeTo(TMatrixDSym(7));
    matrix_ = matrix;
    genRef_ = genRef;
}

const TMatrixDSym & ACParticle::matrix() const { return matrix_; }
const ACGenParticleRef & ACParticle::genRef() const { return genRef_; }

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
const TMatrixDSym ACParticle::vertexError() const {
    TMatrixDSym tmp(3);
    for (int i=0; i!=3; i++) for (int j=0; j!=3; j++) tmp(i,j) = matrix()(i,j);

    return tmp;
}
double ACParticle::svx() const { if (matrix()(0,0)<.0) return -1.0; return sqrt(matrix()(0,0)); }
double ACParticle::svy() const { if (matrix()(1,1)<.0) return -1.0; return sqrt(matrix()(1,1)); }
double ACParticle::svz() const { if (matrix()(2,2)<.0) return -1.0; return sqrt(matrix()(2,2)); }
double ACParticle::spx() const { if (matrix()(3,3)<.0) return -1.0; return sqrt(matrix()(3,3)); }
double ACParticle::spy() const { if (matrix()(4,4)<.0) return -1.0; return sqrt(matrix()(4,4)); }
double ACParticle::spz() const { if (matrix()(5,5)<.0) return -1.0; return sqrt(matrix()(5,5)); }
double ACParticle::sm() const { if (matrix()(6,6)<.0) return -1.0; return sqrt(matrix()(6,6)); }

/// todo: sqoverp, ...
