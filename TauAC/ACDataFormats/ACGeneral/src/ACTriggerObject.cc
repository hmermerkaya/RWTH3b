#include "../interface/ACTriggerObject.h"

ClassImp(ACTriggerObject)

ACTriggerObject::ACTriggerObject() {
    pdgID_ = 0;
    pt_ = 0.;
    eta_ = 0.;
    phi_ = 0.;
    mass_ = 0.;
}
ACTriggerObject::ACTriggerObject(int pdgID, float pt, float eta, float phi, float mass) {
    pdgID_ = pdgID;
    pt_ = pt;
    eta_ = eta;
    phi_ = phi;
    mass_ = mass;
}

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values

float ACTriggerObject::px() const {return pt()*std::cos(phi());}
float ACTriggerObject::py() const {return pt()*std::sin(phi());}
float ACTriggerObject::pz() const {return pt()*std::sinh(eta());}
float ACTriggerObject::p() const {return pt()*std::cosh(eta());}
float ACTriggerObject::energy() const {return std::sqrt(std::pow(mass(),2)+std::pow(p(),2));}
