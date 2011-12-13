#ifndef ACRef_h
#define ACRef_h
/**
 Class to store a reference to an object derived from ACCommon.

 @tparam T Class type of the referenced object (e.g. T=ACParticle).

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCommon.h"

#include "TRef.h"

// forward declarations
class ACCandidate;
class ACGenParticle;
class ACParticle;
class ACFitParticle;
class ACJet;
class ACPFTau;
class ACVertex;
template<class T> class ACDecayBase;
class ACFittedDecay;
class ACFittedThreeProngDecay;

template<class T> class ACRef : public ACCommon {
public:
    virtual std::string classname() const { return "ACRef"; }
    ACRef();
    ACRef(T* obj);
    ACRef(const T* obj);
    virtual ~ACRef() {};

    /// tests whether the reference is valid
    bool isValid() const;
    /// returns C++ pointer to the item
    const T * get() const;
    /// dereference operator
    const T & operator*() const;
    /// member dereference operator
    const T * operator->() const;

    ClassDef(ACRef,1);

protected:
    /// underlying reference object
    TRef ref_;
};

template<class T> ACRef<T>::ACRef() {
    ref_ = TRef();
}
template<class T> ACRef<T>::ACRef(T* obj) {
    ref_ = (TObject*)obj;
}
template<class T> ACRef<T>::ACRef(const T* obj) {
    ref_ = (TObject*)obj;
}
template<class T> bool ACRef<T>::isValid() const {
    return ref_.IsValid();
}
template<class T> const T * ACRef<T>::get() const {
    if (!ref_.IsValid()) std::cerr << "Warning: reference is not valid!" << std::endl;
    return (const T*)(ref_.GetObject());
}
template<class T> const T & ACRef<T>::operator*() const {
    return *get();
}
template<class T> const T * ACRef<T>::operator->() const {
    return get();
}

typedef ACRef<ACCandidate> ACCandidateRef;
typedef ACRef<ACGenParticle> ACGenParticleRef;
typedef ACRef<ACParticle> ACParticleRef;
typedef ACRef<ACFitParticle> ACFitParticleRef;
typedef ACRef<ACJet> ACJetRef;
typedef ACRef<ACPFTau> ACPFTauRef;
typedef ACRef<ACVertex> ACVertexRef;
typedef ACRef<ACDecayBase<ACGenParticle> > ACGenDecayRef;
typedef ACRef<ACDecayBase<ACParticle> > ACDecayRef;
typedef ACRef<ACFittedDecay> ACFittedDecayRef;
typedef ACRef<ACFittedThreeProngDecay> ACFittedThreeProngDecayRef;

#endif
