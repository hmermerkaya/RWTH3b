#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"


const ParticleMass ThreeProngTauCreator::piMass = .13957018;//PYTHIA: 0.140, PDG: 0.13957018
const ParticleMass ThreeProngTauCreator::tauMass = 1.777;

int ThreeProngTauCreator::create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks){
	std::vector<RefCountedKinematicParticle> *daughters = new std::vector<RefCountedKinematicParticle>; //3 particles (3pi)
	std::vector<RefCountedKinematicParticle> *neutrinos = new std::vector<RefCountedKinematicParticle>; //1 or 2 particles due to ambiguity (nuGuess1 + nuGuess2)
	std::vector<reco::TrackRef> input = inputTracks;
	reco::Vertex primVtx = primaryVertex;
	
	if(!createStartScenario(input, *daughters, *neutrinos, primVtx)) return 0;
	if (daughters->size()!=3 ||(neutrinos->size()!=1 && neutrinos->size()!=2)){
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::create: wrong daughter size. found "<<daughters->size()<<" pis and "<<neutrinos->size()<<" nus. Skip this tauCand";
		return 0;
	}
    
//    for ( unsigned int i = 0; i < daughters->size(); ++i ) {
//        std::cerr << "pion_" << i+1 << " matrix: " << std::endl;
//        std::cerr << daughters->at(i)->initialState().kinematicParametersError().matrix() << std::endl;
//    }
//    std::cerr << "neutrino matrix: " << std::endl;
//    std::cerr << neutrinos->at(0)->initialState().kinematicParametersError().matrix() << std::endl;
    
	//in this version createStartScenario always rotates up to thetaMax so that there is always only one solution
	daughters->push_back(neutrinos->at(0));
	
	bool fitWorked = kinematicRefit(*daughters, primVtx);

	delete daughters;
	delete neutrinos;
	
	if(fitWorked) return 1;
	else return 0;
}

bool ThreeProngTauCreator::createStartScenario(std::vector<reco::TrackRef> &input, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos, reco::Vertex &primVtx){
	KinematicParticleFactoryFromTransientTrack kinFactory;
	float piMassSigma = sqrt(pow(10.,-12.));//not to small to avoid singularities
	float piChi = 0., piNdf = 0.;//only initial values
	
	if (input.size()<3) {
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad track size = "<<input.size();
		return false;
	}
	if (input.size()>3){
		if(!choose3bestTracks(input, primVtx)) return false;
	}else{
		if(!sumCharge(input)){
			LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand due to bad charge sum.";
			return false;
		}
		double massA1 = getInvariantMass(input, ThreeProngTauCreator::piMass);
		if(massA1 > 2.0 || massA1 < 3*piMass){//soft upper value
			LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand due to bad a1 mass.";
			return false;
		}
	}
	
	selectedTracks_ = input;
	std::vector<reco::TransientTrack> transTrkVect = convToTransTrck(selectedTracks_);
	if(transTrkVect.size()!=3){
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad transient track size = "<<transTrkVect.size();
		return false;
	}
	TransientVertex secVtx;
	if(!checkSecVtx(transTrkVect, secVtx)){
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Skip tauCand.";
		return false;
	}
	double massA1 = getInvariantMass(selectedTracks_, ThreeProngTauCreator::piMass);
	TLorentzVector lorentzA1 = getSumTLorentzVec(selectedTracks_, massA1);
	VertexRotation vtxC(lorentzA1, 0);
	double thetaMax = vtxC.calcThetaMax();
	if(lorentzA1.M() > 2.0 || lorentzA1.M() < 3*ThreeProngTauCreator::piMass){//soft upper value due to neutrino resolution
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad a1 mass = "<<lorentzA1.M()<<". Skip tauCand.";
		return false;
	}
	
	double theta0;
	TVector3 tauFlghtDir;
	double significance = vtxC.rotatePV(primVtx, secVtx, theta0, tauFlghtDir);//rotation is forced to reach thetaMax
	if(significance > 0.0) modifiedPV_ = primVtx;//if rotation was needed store modified vertex
	
	for(unsigned int i = 0; i!=transTrkVect.size();i++){
		pions.push_back(kinFactory.particle(transTrkVect[i],ThreeProngTauCreator::piMass,piChi,piNdf,secVtx.position(),piMassSigma));
	}
	
	if(theta0>TMath::Pi()/2){
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Unrealistic GJ angle ("<<theta0<<"). tauCand skipped!";
		return false;
	}
	
	if ((thetaMax - theta0) < -pow(10.0,-10.0)*thetaMax && thetaMax>0) {
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: thetaGJ = "<<theta0<<" but replaced by thetaMax = "<<thetaMax;
		theta0 = thetaMax;
	}

	std::pair<double,double> tauSolutions = getTauMomentumMagnitudes(lorentzA1.M(),lorentzA1.P(),ThreeProngTauCreator::tauMass,theta0);//use a pair to prepare for ambiguities
	LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: tau solutions = "<<tauSolutions.first<<", "<<tauSolutions.second;
	bool ambiguity = false;
	if(fabs(tauSolutions.first-tauSolutions.second) > pow(10.0,-6)) ambiguity = true;
	
	TLorentzVector tauGuess1;
	tauFlghtDir = tauFlghtDir.Unit();
	tauGuess1.SetXYZM(tauFlghtDir.X()*tauSolutions.first , tauFlghtDir.Y()*tauSolutions.first , tauFlghtDir.Z()*tauSolutions.first , ThreeProngTauCreator::tauMass);
	neutrinos.push_back(unknownNu(tauGuess1, lorentzA1, secVtx));
	if(ambiguity==true){
		TLorentzVector tauGuess2;
		tauGuess2.SetXYZM(tauFlghtDir.X()*tauSolutions.second, tauFlghtDir.Y()*tauSolutions.second, tauFlghtDir.Z()*tauSolutions.second, ThreeProngTauCreator::tauMass);
		neutrinos.push_back(unknownNu(tauGuess2, lorentzA1, secVtx));
	}
	if(neutrinos.size() != 1 && neutrinos.size() != 2){
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::createStartScenario: Bad neutrino size = "<<neutrinos.size();
		return false;
	}
		
	return true;
}
bool ThreeProngTauCreator::kinematicRefit(std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex &primVtx){
	if(unfitDaughters.size()!=4){
        edm::LogError("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.";
		return false;
	}

	std::vector<MultiTrackKinematicConstraint* > constraintVector;
	MultiTrackKinematicConstraint *tauMass_c = new  MultiTrackMassKinematicConstraint(ThreeProngTauCreator::tauMass, unfitDaughters.size());
	constraintVector.push_back(tauMass_c);
	GlobalPoint linP(primVtx.x(), primVtx.y(), primVtx.z());
//	MultiTrackKinematicConstraint *pointing_c = new MultiTrackPointingKinematicConstraint(linP);
	MultiTrackKinematicConstraint *pointing_c = new MultiTrackVertexLinkKinematicConstraint(linP);
	constraintVector.push_back(pointing_c);
	MultiTrackKinematicConstraint *combiC = new CombinedKinematicConstraint(constraintVector);

	GlobalPoint vtxGuess = unfitDaughters[3]->currentState().globalPosition();//nu was created at common/corrected vertex of pions
	try{
		kinTree_ = kcvFitter_->fit(unfitDaughters, combiC);
	}catch(VertexException){//("KinematicStatePropagator without material::propagation failed!")
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: VertexException. Skip tau candidate.";
		return false;
	}

	delete combiC;
	delete pointing_c;
	delete tauMass_c;
	
	if(kinTree_->isValid()) return true;
	else{
		edm::LogVerbatim("ThreeProngTauCreator")<<"ThreeProngTauCreator::kinematicRefit: ERROR! Tree is not valid. Skip tauCand.";
		return false;
	}
}

bool ThreeProngTauCreator::choose3bestTracks(std::vector<reco::TrackRef> &input, reco::Vertex & pVtx){
	sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);
	std::vector<std::vector<reco::TrackRef> > combis = permuteCombinations(input);
	if(combis.size()==0){
		edm::LogVerbatim("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks: ERROR! Bad combinations. Skip it.";
		return false;
	}
	std::vector<std::pair<int,double> > movements;
	unsigned index=0;
	for (std::vector<std::vector<reco::TrackRef> >::iterator iter=combis.begin(); iter!=combis.end();) {
		if(!sumCharge(*iter)){
			iter = combis.erase(iter);
			LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks: erased combi due to wrong charge sum. "<<combis.size()<<" combis left.";
			continue;
		}
		double massA1 = getInvariantMass(*iter, ThreeProngTauCreator::piMass);
		if(massA1 > 2.0 || massA1 < 3*ThreeProngTauCreator::piMass){//soft upper value
			iter = combis.erase(iter);
			LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks: erased combi due to wrong mass. "<<combis.size()<<" combis left.";
			continue;
		}
		TransientVertex tmpVtx;
		std::vector<reco::TransientTrack> trks = convToTransTrck(*iter);
		if(!checkSecVtx(trks, tmpVtx)){
			iter = combis.erase(iter);
			LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks: erased combi due to bad vertex. "<<combis.size()<<" combis left.";
			continue;
		}
		TLorentzVector lorentzA1 = getSumTLorentzVec(*iter, massA1);
		VertexRotation vtxC(lorentzA1);
		double theta0;
		TVector3 tauFlghtDir;
		reco::Vertex pvTemp = pVtx;//do not modify original pv here
		double significance = vtxC.rotatePV(pvTemp, tmpVtx, theta0, tauFlghtDir);
		movements.push_back(std::make_pair(index,significance));//significance of vertex modification
		
		++index;
		++iter;//only moved if nothing was deleted
	}
	if (combis.size()<1){
		if (combis.size()>=1) LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks:No combi survived.";
		return false;
	}
	
	if (combis.size()>1){//chose the pair with smallest vertex rotation needed
		sort(movements.begin(), movements.end(), pairSecond<int, double>);
		LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::choose3bestTracks:Too much combis ("<<combis.size()<<") left. Take one with smallest vtx correction movement ("<<movements.front().second<<" [sigma]), second best is ("<<movements.at(1).second<<" [sigma]).";
		unsigned int i = movements.front().first;
		input.assign( combis.at(i).begin(), combis.at(i).end() );
		
	}else input.assign( combis.front().begin(), combis.front().end() );
	
	return true;
}
bool ThreeProngTauCreator::sumCharge(std::vector<reco::TrackRef> &input){
	int sum = abs(input.at(0)->charge() + input.at(1)->charge() + input.at(2)->charge());
	if(sum == 1) return true;
	return false;
}
std::vector<reco::TransientTrack> ThreeProngTauCreator::convToTransTrck(std::vector<reco::TrackRef> &input){
	//TransientTrackBuilder delivers reco::TransientTrack*
	std::vector<reco::TransientTrack> transTrkVct;
	for (std::vector<reco::TrackRef>::iterator iter=input.begin(); iter!=input.end(); ++iter) {
		transTrkVct.push_back( transTrackBuilder_.build( *iter ) );
	}
	return transTrkVct;
}
bool ThreeProngTauCreator::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx){
	if(trkVct.size()<2){
		LogTrace("ThreeProngTauCreator")<<"Can't check SecVertex: Only "<<trkVct.size()<<" Tracks.";
		return false;
	}else{
		bool useAdaptive = false;
		if(useAdaptive){
			AdaptiveVertexFitter avf;
			avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
			try{
				transVtx = avf.vertex(trkVct); //AdaptiveVertexFitter
			}catch(...){
				LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}else{
			KalmanVertexFitter kvf(true);
			try{
				transVtx = kvf.vertex(trkVct); //KalmanVertexFitter
			}catch(...){
				LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex fit failed. Skip it.";
				return false;
			}
		}
		if(!transVtx.isValid()) LogTrace("ThreeProngTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary vertex not valid.";
		if(!useAdaptive){
			if(!transVtx.hasRefittedTracks()){
				LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has 0 refitted tracks.";
				return false;
			}else if(transVtx.refittedTracks().size()!=trkVct.size()){
				LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::checkSecVtx: Secondary has only "<<transVtx.refittedTracks().size()<<" refitted of "<<trkVct.size()<<" initial tracks.";
				return false;
			}
		}
		
		return transVtx.isValid();
	}
}
std::pair<double,double> ThreeProngTauCreator::getTauMomentumMagnitudes(double ma1,double pa1,double M,double theta){
	double root = sqrt((pow(ma1,2) + pow(pa1,2))*(pow(pow(ma1,2) - pow(M,2),2) -4*pow(M,2)*pow(pa1,2)*pow(sin(theta),2))),
	numerator = (pow(ma1,2) + pow(M,2))*pa1*cos(theta),
	denominator = (2*pow(ma1,2) + 2*pow(pa1,2)*pow(sin(theta),2));
	double tauMomentumMagnitude1 = (numerator - root) / denominator,
	tauMomentumMagnitude2 = (numerator + root)	/ denominator;
	
	if(!(tauMomentumMagnitude1>0)) tauMomentumMagnitude1 = numerator / denominator; //catch 'nan' and replace by pseudo mean value
	if(!(tauMomentumMagnitude2>0)) tauMomentumMagnitude2 = tauMomentumMagnitude1;
	
	if(!(tauMomentumMagnitude1>0)) tauMomentumMagnitude1 = 0; //catch 'nan'
	if(!(tauMomentumMagnitude2>0)) tauMomentumMagnitude2 = 0;
	
	return std::make_pair(tauMomentumMagnitude1,tauMomentumMagnitude2);
}
RefCountedKinematicParticle ThreeProngTauCreator::unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx){
	GlobalVector nuImpGuess(1,1,1);
	
	TLorentzVector nuGuess = tauGuess-a1;
	LogTrace("KinematicTauCreator")<<"ThreeProngTauCreator::unknownNu: nuGuess (vx, vy, vz, px,py,pz,m) "<<secVtx.position().x()<<","<<secVtx.position().y()<<","<<secVtx.position().z()<<","<<nuGuess.Px()<<","<<nuGuess.Py()<<","<<nuGuess.Pz()<<","<<nuGuess.M();
	if(tauGuess.P()==0)  nuGuess.SetXYZM(1,1,1,0);
	nuImpGuess = GlobalVector(nuGuess.X(),nuGuess.Y(),nuGuess.Z());
	return virtualKinematicParticle(secVtx, nuImpGuess);
}
RefCountedKinematicParticle ThreeProngTauCreator::virtualKinematicParticle(TransientVertex & vtxGuess, GlobalVector impulsGuess){
	VirtualKinematicParticleFactory factory;
	//(x,y,z,p_x,p_y,p_z,m)
	const KinematicParameters parameters(AlgebraicVector7(vtxGuess.position().x(),vtxGuess.position().y(),vtxGuess.position().z(),impulsGuess.x(),impulsGuess.y(),impulsGuess.z(),pow(10.,-10.)));//Use MET as initial guess for pt?
	ROOT::Math::SVector<double,28> svector28;
	for(unsigned int i=1; i!=22; i++) svector28(i-1) = 0.0;
	for(unsigned int i=22; i!=28; i++) svector28(i-1) = 0.0;//correlation between mass and momentum/vertex
	for(unsigned int n=1; n!=7; n++) svector28(n*(n+1)/2 - 1) = pow(10.,2.);//diagonals, huge error method
	svector28(27) = pow(10.,-12.);//mass error
    //insert 3prong vertex errors
    svector28[ 0] = vtxGuess.positionError().cxx();
    svector28[ 1] = vtxGuess.positionError().cyx();
    svector28[ 2] = vtxGuess.positionError().cyy();
    svector28[ 3] = vtxGuess.positionError().czx();
    svector28[ 4] = vtxGuess.positionError().czy();
    svector28[ 5] = vtxGuess.positionError().czz();
    
    if (std::abs(impulsGuess.x()) >= 5.0) {
        svector28[ 9] = pow(impulsGuess.x(), 2); //assume an error of 100% of the neutrino momentum component
    } else {
        svector28[ 9] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
    }
    if (std::abs(impulsGuess.y()) >= 5.0) {
        svector28[14] = pow(impulsGuess.y(), 2); //assume an error of 100% of the neutrino momentum component
    } else {
        svector28[14] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
    }
    if (std::abs(impulsGuess.z()) >= 5.0) {
        svector28[20] = pow(impulsGuess.z(), 2); //assume an error of 100% of the neutrino momentum component
    } else {
        svector28[20] = pow(5.0, 2); //for momenta smaller than 5 GeV set the error to a static value
    }
	
	ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > matrix(svector28);
	
	const KinematicParametersError parametersError(matrix);
	const TrackCharge charge = 0;
	KinematicState kineState(parameters, parametersError, charge, transTrackBuilder_.field());
	float chiSquared=0.0, degreesOfFr=0.0;
	
	return factory.particle(kineState, chiSquared, degreesOfFr, 0,0);
}
template <typename T> std::vector<std::vector<T> > ThreeProngTauCreator::permuteCombinations(const std::vector<T> &vect){
	std::vector<std::vector<T> > combis;
	typename std::vector<T>::const_iterator iter1, iter2, iter3;
	for (iter1 = vect.begin(); iter1!=vect.end(); ++iter1) {
		iter2 = iter1;
		++iter2;
		for (; iter2!=vect.end(); ++iter2) {
			iter3 = iter2;
			++iter3;
			for (; iter3!=vect.end(); ++iter3) {
				std::vector<T> newCombi;
				newCombi.push_back(*iter1);
				newCombi.push_back(*iter2);
				newCombi.push_back(*iter3);
				combis.push_back(newCombi);
			}			
		}
	}
	return combis;
}
