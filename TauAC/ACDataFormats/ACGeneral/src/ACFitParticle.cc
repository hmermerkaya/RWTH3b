#include "../interface/ACFitParticle.h"

ClassImp(ACFitParticle)

ACFitParticle::ACFitParticle():
ACParticle() {
    initialP4_.SetXYZM(-1000,-1000,-1000,-1000);
    initialVertex_.SetXYZ(-1000,-1000,-1000);
    initialMatrix_.ResizeTo(TMatrixDSym(7));
}
ACFitParticle::ACFitParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex, const TMatrixDSym & matrix, const TLorentzVector & initialP4, const TVector3 & initialVertex, const TMatrixDSym & initialMatrix, ACGenParticleRef genRef, ACGenDecayRef genDecayRef):
ACParticle(status, name, pdgId, charge, p4, vertex, matrix, genRef, genDecayRef) {
    initialP4_ = initialP4;
    initialVertex_ = initialVertex;
    initialMatrix_.ResizeTo(TMatrixDSym(7));
    initialMatrix_ = initialMatrix;
}

const TLorentzVector & ACFitParticle::initialP4() const { return initialP4_; }
const TVector3 & ACFitParticle::initialVertex() const { return initialVertex_; }
const TMatrixDSym & ACFitParticle::initialMatrix() const { return initialMatrix_; }

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
const ACVertex ACFitParticle::initialVtx() const {
    return ACVertex(initialVertex(), true, initialVertexError());
}
double ACFitParticle::pInitial() const { return initialP4().P(); }
double ACFitParticle::pxInitial() const { return initialP4().Px(); }
double ACFitParticle::pyInitial() const { return initialP4().Py(); }
double ACFitParticle::pzInitial() const { return initialP4().Pz(); }
double ACFitParticle::ptInitial() const { return initialP4().Pt(); }

const TMatrixDSym ACFitParticle::initialVertexError() const {
    TMatrixDSym tmp(3);
    for (int i=0; i!=3; i++) for (int j=0; j!=3; j++) tmp(i,j) = initialMatrix()(i,j);

    return tmp;
}
double ACFitParticle::svxInitial() const { if (initialMatrix()(0,0)<.0) return -1.0; return sqrt(initialMatrix()(0,0)); }
double ACFitParticle::svyInitial() const { if (initialMatrix()(1,1)<.0) return -1.0; return sqrt(initialMatrix()(1,1)); }
double ACFitParticle::svzInitial() const { if (initialMatrix()(2,2)<.0) return -1.0; return sqrt(initialMatrix()(2,2)); }
double ACFitParticle::spxInitial() const { if (initialMatrix()(3,3)<.0) return -1.0; return sqrt(initialMatrix()(3,3)); }
double ACFitParticle::spyInitial() const { if (initialMatrix()(4,4)<.0) return -1.0; return sqrt(initialMatrix()(4,4)); }
double ACFitParticle::spzInitial() const { if (initialMatrix()(5,5)<.0) return -1.0; return sqrt(initialMatrix()(5,5)); }
double ACFitParticle::smInitial() const { if (initialMatrix()(6,6)<.0) return -1.0; return sqrt(initialMatrix()(6,6)); }

double ACFitParticle::sptInitial() const {
    if(spxInitial()<.0 || spyInitial()<.0) return -1.0;
    return sqrt(pow(pxInitial()*spxInitial(),2)+pow(pyInitial()*spyInitial(),2))/ptInitial();
}
double ACFitParticle::spInitial() const {
    if(spxInitial()<.0 || spyInitial()<.0 || spzInitial()<.0) return -1.0;
    return sqrt(pow(pxInitial()*spxInitial(),2)+pow(pyInitial()*spyInitial(),2)+pow(pzInitial()*spzInitial(),2))/pInitial();
}
double ACFitParticle::qoverpInitial() const {
    return charge()/initialP4().P();
}
double ACFitParticle::lambdaInitial() const {
    return atan(initialP4().Pz()/initialP4().Pt());
}
double ACFitParticle::dxyInitial() const {
    return -initialVertex().X()*sin(initialP4().Phi())+initialVertex().Y()*cos(initialP4().Phi());
}
double ACFitParticle::dszInitial() const {
    return initialVertex().Z()*cos(lambda())-(initialVertex().X()*cos(initialP4().Phi())+initialVertex().Y()*sin(initialP4().Phi()))*sin(lambda());
}
