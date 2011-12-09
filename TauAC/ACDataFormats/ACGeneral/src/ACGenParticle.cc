#include "../interface/ACGenParticle.h"

ClassImp(ACGenParticle)

ACGenParticle::ACGenParticle():
ACCandidate() {
    status_ = -1;
    name_ = "";
    pdgId_ = 0;
    charge_ = -1000.;
    vertex_.SetXYZ(-1000,-1000,-1000);
    genDecayRef_ = ACGenDecayRef();
}
ACGenParticle::ACGenParticle(int status, std::string name, int pdgId, float charge, const TLorentzVector & p4, const TVector3 & vertex, const ACGenDecayRef & genDecayRef):
ACCandidate(p4, charge) {
    status_ = status;
    name_ = name;
    pdgId_ = pdgId;
    vertex_ = vertex;
    genDecayRef_ = genDecayRef;
}

int ACGenParticle::status() const { return status_; }
std::string ACGenParticle::name() const { return name_; }
int ACGenParticle::pdgId() const { return pdgId_; }
const TVector3 & ACGenParticle::vertex() const { return vertex_; }
const ACGenDecayRef & ACGenParticle::genDecayRef() const { return genDecayRef_; }

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
double ACGenParticle::vx() const { return vertex().X(); }
double ACGenParticle::vy() const { return vertex().Y(); }
double ACGenParticle::vz() const { return vertex().Z(); }

double ACGenParticle::qoverp() const { return charge()/p(); }
double ACGenParticle::lambda() const { return atan(pz()/pt()); }
double ACGenParticle::dxy() const { return -vx()*sin(phi())+vy()*cos(phi()); }
double ACGenParticle::dsz() const { return vz()*cos(lambda())-(vx()*cos(phi())+vy()*sin(phi()))*sin(lambda()); }
