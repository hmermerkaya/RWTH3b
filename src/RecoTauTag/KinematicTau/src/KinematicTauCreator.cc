#include "RecoTauTag/KinematicTau/interface/KinematicTauCreator.h"


KinematicTauCreator::KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder):
transTrackBuilder_(transTrackBuilder)
{
    kcvFitter_ = new KinematicConstrainedVertexFitter();
    edm::ParameterSet defaultConfig;
    defaultConfig.addParameter("maxDistance", .001);
    defaultConfig.addParameter("maxNbrOfIterations", 20);
    defaultConfig.addParameter("maxOfInitialValue", 9999.);
    kcvFitter_->setParameters(defaultConfig);
}

KinematicTauCreator::KinematicTauCreator(const TransientTrackBuilder & transTrackBuilder, const edm::ParameterSet& cfg):
transTrackBuilder_(transTrackBuilder)
{
    kcvFitter_ = new KinematicConstrainedVertexFitter();
    kcvFitter_->setParameters(cfg);
}

KinematicTauCreator::~KinematicTauCreator()
{
    delete kcvFitter_;
}

reco::PFTau KinematicTauCreator::getPFTau()
{
    math::XYZTLorentzVector ptau;
    for ( unsigned int i = 0; i < kinTree_->daughterParticles().size(); i++ ) {
        if (std::abs(int((kinTree_->daughterParticles().at(i))->currentState().particleCharge())) == 1) {
            ptau += math::XYZTLorentzVector((kinTree_->daughterParticles().at(i))->currentState().globalMomentum().x(),
                                            (kinTree_->daughterParticles().at(i))->currentState().globalMomentum().y(),
                                            (kinTree_->daughterParticles().at(i))->currentState().globalMomentum().z(),
                                            sqrt((kinTree_->daughterParticles().at(i))->currentState().globalMomentum().mag2()+(kinTree_->daughterParticles().at(i))->currentState().mass()*(kinTree_->daughterParticles().at(i))->currentState().mass()));
        }
    }
    
    return reco::PFTau(int((kinTree_->topParticle()->currentState()).particleCharge()), ptau, (math::XYZPoint)((kinTree_->topParticle()->currentState()).globalPosition()));
}

std::vector<math::XYZTLorentzVector> KinematicTauCreator::getRefittedChargedHadrons()
{
    std::vector<math::XYZTLorentzVector> tmpvec;
    for ( unsigned int i = 0; i < kinTree_->daughterParticles().size(); i++ ) {
        if (std::abs(int((kinTree_->daughterParticles().at(i))->currentState().particleCharge())) == 1) {
            tmpvec.push_back(math::XYZTLorentzVector((kinTree_->daughterParticles().at(i))->currentState().globalMomentum().x(),
                                            (kinTree_->daughterParticles().at(i))->currentState().globalMomentum().y(),
                                            (kinTree_->daughterParticles().at(i))->currentState().globalMomentum().z(),
                                            sqrt((kinTree_->daughterParticles().at(i))->currentState().globalMomentum().mag2()+(kinTree_->daughterParticles().at(i))->currentState().mass()*(kinTree_->daughterParticles().at(i))->currentState().mass())));
        }
    }
    
    return tmpvec;
}

std::vector<reco::TrackRef> KinematicTauCreator::getSelectedTracks()
{
    return selectedTracks_;
}


RefCountedKinematicTree KinematicTauCreator::getKinematicTree()
{
    return kinTree_;
}
