#ifndef ACDecayBase_h
#define ACDecayBase_h
/**
 Base class to store decay information. The template argument expects either T=ACParticle or T=ACFitParticle

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "ACCommon.h"
#include "ACParticle.h"
#include "ACFitParticle.h"
#include "ACRef.h"

template<class T> class ACDecayBase : public ACCommon {
public:
    virtual std::string classname() const { return "ACDecayBase"; }
    ACDecayBase();
    ACDecayBase(std::vector<ACRef<T> > particles);
    ACDecayBase(std::vector<ACRef<T> > particles, std::string name, ACGenDecayRef genDecayRef = ACGenDecayRef());
    virtual ~ACDecayBase() {}

    /// all particles
    const std::vector<ACRef<T> > * particles() const;
    const std::string name() const;
    const ACGenDecayRef & genDecayRef() const;

    /// the first entry of the particle collection (assuming mother particle)
    const T* topParticle() const;
    /// all daughters (ignoring the top particle)
    void daughters(std::vector<const T*>& daughters) const;
    /// all charged daughters (ignoring the top particle and neutrals)
    void chargedDaughters(std::vector<const T*>& chargedDaughters) const;
    /// all neutral daughters (ignoring the top particle and charged ones)
    void neutralDaughters(std::vector<const T*>& neutralDaughters) const;
    /// all pi0 involved in the decay
    void pi0Daughters(std::vector<const T*>& pi0Daughters) const;
    /// number of particles assigned to this decay
    int numberOfParticles() const;
    /// number of daughters (ignoring the top particle)
    int numberOfDaughters() const;
    /// number of charged daughters (ignoring the top particle and neutrals)
    int numberOfChargedDaughters() const;
    /// number of neutral daughters (ignoring the top particle and charged ones)
    int numberOfNeutralDaughters() const;
    /// number of pi0s
    int numberOfPi0Daughters() const;
    /// check for four vector conservation
    double conservation() const;
    /// assumes a tau decay and returns the decay mode according to the TauPOG definition: mode = 5 * (chargedDaughters - 1) + number of pi0's
    int decayMode();

    ClassDef(ACDecayBase,1);

protected:
    /// particle collection (first one is assumed to correspond to the mother particle)
    std::vector<ACRef<T> > particles_;
    /// decay name
    std::string name_;
    /// reference to the generated decay
    ACGenDecayRef genDecayRef_;
};

typedef ACDecayBase<ACParticle> ACDecay; // this is just a typedef until we need a class with additional functionality
typedef ACDecay* pACDecay;

template<class T> ACDecayBase<T>::ACDecayBase() {
    particles_.clear();
    name_ = "";
    genDecayRef_ = ACGenDecayRef();
}
template<class T> ACDecayBase<T>::ACDecayBase(std::vector<ACRef<T> > particles) {
    particles_ = particles;
    name_ = "";
    genDecayRef_ = ACGenDecayRef();
}
template<class T> ACDecayBase<T>::ACDecayBase(std::vector<ACRef<T> > particles, std::string name, ACGenDecayRef genDecayRef) {
    particles_ = particles;
    name_ = name;
    genDecayRef_ = genDecayRef;
}

template<class T> const std::vector<ACRef<T> > * ACDecayBase<T>::particles() const {
    return &particles_;
}
template<class T> const std::string ACDecayBase<T>::name() const { return name_; }
template<class T> const ACGenDecayRef & ACDecayBase<T>::genDecayRef() const { return genDecayRef_; }

template<class T> const T* ACDecayBase<T>::topParticle() const {
    if (particles()->size() < 1) return 0;
    if (!particles()->front().isValid()) return 0;
    return particles()->front().get();
}

template<class T> void ACDecayBase<T>::daughters(std::vector<const T*>& daughters) const {
    daughters.clear();
    for (typename std::vector<ACRef<T> >::const_iterator iter = particles()->begin(); iter != particles()->end(); ++iter) {
        if (iter == particles()->begin()) continue;
        daughters.push_back(iter->get());
    }
}
template<class T> void ACDecayBase<T>::chargedDaughters(std::vector<const T*>& chargedDaughters) const {
    chargedDaughters.clear();
    for (typename std::vector<ACRef<T> >::const_iterator iter = particles()->begin(); iter != particles()->end(); ++iter) {
        if (iter == particles()->begin() || (iter->get())->charge() == 0.0) continue;
        chargedDaughters.push_back(iter->get());
    }
}
template<class T> void ACDecayBase<T>::neutralDaughters(std::vector<const T*>& neutralDaughters) const {
    neutralDaughters.clear();
    for (typename std::vector<ACRef<T> >::const_iterator iter = particles()->begin(); iter != particles()->end(); ++iter) {
        if (iter == particles()->begin() || (iter->get())->charge() != 0.0) continue;
        neutralDaughters.push_back(iter->get());
    }
}
template<class T> void ACDecayBase<T>::pi0Daughters(std::vector<const T*>& pi0Daughters) const {
    pi0Daughters.clear();
    for (typename std::vector<ACRef<T> >::const_iterator iter = particles()->begin(); iter != particles()->end(); ++iter) {
        if (iter == particles()->begin() || abs((iter->get())->pdgId()) != 111) continue;
        pi0Daughters.push_back(iter->get());
    }
}
template<class T> int ACDecayBase<T>::numberOfParticles() const { return particles()->size(); }
template<class T> int ACDecayBase<T>::numberOfDaughters() const {
    std::vector<const T*> tmp;
    daughters(tmp);
    return tmp.size();
}
template<class T> int ACDecayBase<T>::numberOfChargedDaughters() const {
    std::vector<const T*> tmp;
    chargedDaughters(tmp);
    return tmp.size();
}
template<class T> int ACDecayBase<T>::numberOfNeutralDaughters() const {
    std::vector<const T*> tmp;
    neutralDaughters(tmp);
    return tmp.size();
}
template<class T> int ACDecayBase<T>::numberOfPi0Daughters() const {
    std::vector<const T*> tmp;
    pi0Daughters(tmp);
    return tmp.size();
}
template<class T> double ACDecayBase<T>::conservation() const {
    TLorentzVector mother, sum;
    if (numberOfParticles()<2) {
        printf("Decay::conservation: Error. Bad particle size of %d.\n", numberOfParticles());
        return -1000.0;
    }
    mother = topParticle()->p4();
    for (typename std::vector<ACRef<T> >::const_iterator p = particles()->begin(); p != particles()->end(); ++p) {
        if (p==particles()->begin()) continue;
        sum += (*p)->p4();
    }

    return sum.M() - mother.M();
}
template<class T> int ACDecayBase<T>::decayMode() {
    // decay mode = 5*(number of charged hadrons-1) + number of pi0's.
    int mode = 5 * (numberOfChargedDaughters() - 1) + numberOfPi0Daughters();

    return mode;
}

#endif
