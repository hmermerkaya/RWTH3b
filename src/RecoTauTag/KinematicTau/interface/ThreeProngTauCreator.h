#include "RecoTauTag/KinematicTau/KinematicTauCreator.h"



class ThreeProngTauCreator : public KinematicTauCreator
{
public:
    
private:
    int create(const reco::Vertex& primvtx_, const std::vector<reco::TrackRef>& inputTracks);

    
};