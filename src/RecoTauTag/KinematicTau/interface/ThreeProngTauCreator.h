#include "RecoTauTag/KinematicTau/interface/KinematicTauCreator.h"
#include "RecoTauTag/KinematicTau/interface/VertexRotation.h"

#include <TLorentzVector.h>
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include <RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h>
//kinematic fit:
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include <RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h>
#include <RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h>
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
//own KinematicFit classes
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"


class ThreeProngTauCreator : public KinematicTauCreator
{
public:
	ThreeProngTauCreator(const TransientTrackBuilder & transTrackBuilder):KinematicTauCreator(transTrackBuilder){}
    ThreeProngTauCreator(const TransientTrackBuilder & transTrackBuilder, const edm::ParameterSet& cfg):KinematicTauCreator(transTrackBuilder, cfg){}

private:
    virtual int create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks);

	bool createStartScenario(std::vector<reco::TrackRef> &input, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos, reco::Vertex &primVtx);
	bool kinematicRefit(std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex &primVtx);
	bool choose3bestTracks(std::vector<reco::TrackRef> &input, reco::Vertex & pVtx);
	bool sumCharge(std::vector<reco::TrackRef> &input);
	std::vector<reco::TransientTrack> convToTransTrck(std::vector<reco::TrackRef> &input);
	bool checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx);
	GlobalVector calcPVSVDir(reco::Vertex &primVtx, TransientVertex & transVtx);
	std::pair<double,double> getTauMomentumMagnitudes(double ma1,double pa1,double M,double theta);
	RefCountedKinematicParticle virtualKinematicParticle(GlobalPoint vtxGuess, GlobalVector impulsGuess);	
	RefCountedKinematicParticle unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx);
	template <typename T> std::vector<std::vector<T> > permuteCombinations(const std::vector<T> &vect);
	
	edm::InputTag primVtx_, usedTauCandidatesTag_, inputCollectionTag_;
	int verbosity_;
	unsigned int cnt, cntFound;

	
	int factorial(int num){
		if(num > 1)	return num*factorial(num-1);
		return 1;
	}
	int nOverM(const int& n, const int& m){
		return factorial(n)/(factorial(m)*factorial(n-m));
	}
	template <typename T> double getInvariantMass(const T& tracks, const double mass = 0.140){//if second argument empty default pion is supposed
		double SumPx = 0;
		double SumPy = 0;
		double SumPz = 0;
		double SumE = 0;
		
		for(unsigned int i=0; i<tracks.size(); i++){
			SumPx += tracks[i]->px();
			SumPy += tracks[i]->py();
			SumPz += tracks[i]->pz();
			SumE += sqrt(pow(tracks[i]->p(),2)+pow(mass,2));
		}
		double invmass = sqrt(pow(SumE,2)-pow(SumPx,2)-pow(SumPy,2)-pow(SumPz,2));
		return invmass;
	}
	template <class T> TLorentzVector getSumTLorentzVec(const T& tracks, const double massConstraint){
		double sumPx=0, sumPy=0, sumPz=0;
		for(unsigned int i=0; i!= tracks.size(); i++){
			sumPx += tracks.at(i)->px();
			sumPy += tracks.at(i)->py();
			sumPz += tracks.at(i)->pz();
		}
		TLorentzVector lorentz;
		lorentz.SetXYZM(sumPx, sumPy, sumPz, massConstraint);
		
		return lorentz;
	}
	template <typename S, typename T> static bool pairSecond(const std::pair<S, T> &a, const std::pair<S, T> &b){
		return a.second < b.second;
	}
	template <typename T> static bool cmpPt(const T &a, const T &b){
		return a->pt() > b->pt();
	}
	template <typename T> static bool cmpChi2(const T &a, const T &b){
		return a->chi2() > b->chi2();
	}
};
