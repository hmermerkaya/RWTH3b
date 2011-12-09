#include "../interface/ACFittedThreeProngDecay.h"

ClassImp(ACFittedThreeProngDecay)

ACFittedThreeProngDecay::ACFittedThreeProngDecay():ACFittedDecay() {
    constraints_ = -1;
}
ACFittedThreeProngDecay::ACFittedThreeProngDecay(std::vector<ACFitParticleRef> particles, std::string name, int iterations, float csum, int maxiterations, float mincsum, float chi2, ACPFTauRef PFTauRef, int constraints):
ACFittedDecay(particles, name, iterations, csum, maxiterations, mincsum, chi2, PFTauRef) {
    constraints_ = constraints;
}

int ACFittedThreeProngDecay::ndf() const {
    if (!testDecay()) return -1000;
    int freeParameters = 28-15-4-3;
    int ndf = constraints_ - freeParameters;
    if (ndf != 2) printf("ACFittedThreeProngDecay::ndf: Warning! Unexpected ndf of %i. Expected 2.\n", ndf);

    return ndf;
}
float ACFittedThreeProngDecay::normalizedChi2() const {
    if (!testDecay()) return -1.;
    return (ndf()==0) ? -1 : (float)chi2() / (float)ndf();
}

const ACFitParticle * ACFittedThreeProngDecay::tau() const { if (testDecay()) return particles()->at(0).get(); return 0; }
const ACFitParticle * ACFittedThreeProngDecay::prong0() const { if (testDecay()) return particles()->at(1).get(); return 0; }
const ACFitParticle * ACFittedThreeProngDecay::prong1() const { if (testDecay()) return particles()->at(2).get(); return 0; }
const ACFitParticle * ACFittedThreeProngDecay::prong2() const { if (testDecay()) return particles()->at(3).get(); return 0; }
const ACFitParticle * ACFittedThreeProngDecay::nu() const { if (testDecay()) return particles()->at(4).get(); return 0; }

TLorentzVector ACFittedThreeProngDecay::a1Lorentz(bool prefit) const {
    if (!testDecay()) return TLorentzVector(0,0,0,0);
    TLorentzVector sum;
    for (std::vector<ACFitParticleRef>::const_iterator p = particles()->begin(); p != particles()->end(); ++p) {
        if (p == particles()->begin()) continue; // ignore tau
        if ((*p)->charge()==0.) continue; // ignore nu
        if (!prefit) { // default
            sum += (*p)->p4();
        } else {
            sum += (*p)->initialP4();
        }
    }

    return sum;
}
double ACFittedThreeProngDecay::gjAngle(bool prefit) const {
    if (!testDecay()) return -1.;
    return a1Lorentz(prefit).Angle(tau()->p4().Vect());
}
double ACFittedThreeProngDecay::maxGJAngle(bool prefit) const {
    if (!testDecay()) return -1.;
    TLorentzVector a1 = a1Lorentz(prefit);
    double Mtau = 1.777;
    double argument = (-pow(a1.M(),2.) + pow(Mtau,2.))/(2.*Mtau*a1.P());
    if (fabs(argument) >  1.0) {
//        printf("ACFittedThreeProngDecay::maxGJAngle: Warning! arcsin(%f) = %f. (pa1 %f, ma1 %f)\n", argument, asin(argument), a1.P(), a1.M());
        if (argument >  1.0) argument =  1.0;
        if (argument < -1.0) argument = -1.0;
    }

    return asin(argument);
}

bool ACFittedThreeProngDecay::testDecay() const {
    if (numberOfParticles()!=5) {
        printf("ACFittedThreeProngDecay::testDecay: Bad decay assumption of 5 particles. Found instead %d particles.\n", numberOfParticles());
        return false;
    }
    return true;
}
