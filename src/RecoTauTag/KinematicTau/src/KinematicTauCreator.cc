#include "RecoTauTag/KinematicTau/KinematicTauCreator.h"


KinematicTauCreator::KinematicTauCreator()
{
    kcvFitter = KinematicConstrainedVertexFitter();
}

KinematicTauCreator::KinematicTauCreator(const edm::ParameterSet& cfg)
{
    kcvFitter = KinematicConstrainedVertexFitter();
    kcvFitter.setParameters(cfg);
}

KinematicTauCreator::~KinematicTauCreator()
{
    
}

reco::PFTau KinematicTauCreator::getPFTau()
{
    math::XYZTLorentzVector ptau;
    for ( unsigned int i = 0; i < kinTree->daughterParticles().size(); i++ ) {
        if (std::abs(int((kinTree->daughterParticles().at(i))->currentState().particleCharge())) == 1) {
            ptau += math::XYZTLorentzVector((kinTree->daughterParticles().at(i))->currentState().globalMomentum().x(),
                                            (kinTree->daughterParticles().at(i))->currentState().globalMomentum().y(),
                                            (kinTree->daughterParticles().at(i))->currentState().globalMomentum().z(),
                                            sqrt((kinTree->daughterParticles().at(i))->currentState().globalMomentum().mag2()+(kinTree->daughterParticles().at(i))->currentState().globalMomentum().mass()*(kinTree->daughterParticles().at(i))->currentState().globalMomentum().mass()));
        }
    }
    
    return PFTau(int((kinTree->topParticle()->currentState()).particleCharge()), ptau, (kinTree->topParticle()->currentState()).globalPosition());
}

std::vector<math::XYZTLorentzVector> KinematicTauCreator::getRefittedChargedHadrons()
{
    std::vector<math::XYZTLorentzVector> tmpvec;
    for ( unsigned int i = 0; i < kinTree->daughterParticles().size(); i++ ) {
        if (std::abs(int((kinTree->daughterParticles().at(i))->currentState().particleCharge())) == 1) {
            tmpvec.push_back(math::XYZTLorentzVector((kinTree->daughterParticles().at(i))->currentState().globalMomentum().x(),
                                            (kinTree->daughterParticles().at(i))->currentState().globalMomentum().y(),
                                            (kinTree->daughterParticles().at(i))->currentState().globalMomentum().z(),
                                            sqrt((kinTree->daughterParticles().at(i))->currentState().globalMomentum().mag2()+(kinTree->daughterParticles().at(i))->currentState().globalMomentum().mass()*(kinTree->daughterParticles().at(i))->currentState().globalMomentum().mass())));
        }
    }
    
    return tmpvec;
}

RefCountedKinematicTree KinematicTauCreator::getKinematicTree()
{
    return kinTree;
}

int KinematicTauCreator::iterations()
{
    return iterations;
}

float KinematicTauCreator::csum()
{
    return csum;
}
