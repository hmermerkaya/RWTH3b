#ifndef ACTriggerObject_h
#define ACTriggerObject_h
/**
 Class to store basic properties for trigger objects (according to CMSSW's definition,  optimized in storage size).

 @author Lars Perchalla
 @date 2012
 */

#include <TLorentzVector.h>

#include "ACCommon.h"

class ACTriggerObject : public ACCommon {
public:
    virtual std::string classname() const { return "ACTriggerObject"; }
    ACTriggerObject();
    ACTriggerObject(int pdgID, float pt, float eta, float phi, float mass);
    virtual ~ACTriggerObject() {};


    int pdgID() const { return pdgID_; }
    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float mass() const { return mass_; }

    float px() const ;
    float py() const ;
    float pz() const ;
    float p() const ;
    float energy() const ;
    
    ClassDef(ACTriggerObject,1);

protected:
    /// PDG ID
    int pdgID_;
    /// transversal momentum
    float pt_;
    /// eta
    float eta_;
    /// phi
    float phi_;
    /// mass
    float mass_;    
};

typedef ACTriggerObject* pACTriggerObject;

#endif
