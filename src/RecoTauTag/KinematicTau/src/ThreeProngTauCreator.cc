#include "RecoTauTag/KinematicTau/interface/ThreeProngTauCreator.h"


int ThreeProngTauCreator::create(const reco::Vertex& primaryVertex, const std::vector<reco::TrackRef>& inputTracks){
	std::vector<RefCountedKinematicParticle> *pions = new std::vector<RefCountedKinematicParticle>;//3 particles (3pi)
	std::vector<RefCountedKinematicParticle> *neutrinos = new std::vector<RefCountedKinematicParticle>;//1 or 2 particles due to ambiguity (nuGuess1 + nuGuess2)
	std::vector<reco::TrackRef> input = inputTracks;
	reco::Vertex primVtx = primaryVertex;
	
	if(!createStartScenario(input, *pions, *neutrinos, primVtx)) return 0;
	if (pions->size()!=3 ||(neutrinos->size()!=1 && neutrinos->size()!=2)){
		if(verbosity_>=1) printf("ThreeProngTauCreator::create: wrong daughter size. found %d pis and %d nus. Skip this tauCand\n", pions->size(), neutrinos->size());
		return 0;
	}
	//in this version createStartScenario always rotates up to thetaMax so that there is always only one solution
	pions->push_back(neutrinos->at(0));

	bool fitWorked = kinematicRefit(*pions, primVtx);
	delete pions;
	delete neutrinos;
	
	if(fitWorked) return 1;
	else return 0;
}

bool ThreeProngTauCreator::createStartScenario(std::vector<reco::TrackRef> &input, std::vector<RefCountedKinematicParticle> &pions, std::vector<RefCountedKinematicParticle> &neutrinos, reco::Vertex &primVtx){
	KinematicParticleFactoryFromTransientTrack kinFactory;
	ParticleMass piMass = .140;//PYTHIA: 0.140, PDG: 0.13957018
	float piMassSigma = sqrt(pow(10.,-12.));//not to small to avoid singularities
	float piChi = 0., piNdf = 0.;//only initial values
	
	if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: pion size = %d.\n", input.size());
	if (input.size()<3) {
		if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: Bad track size = %d.\n", input.size());
		return false;
	}
	if (input.size()>3)	if(!choose3bestTracks(input, primVtx)) return false;
	else if(!sumCharge(input)){
		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: Skip tauCand due to bad charge sum.\n");
		return false;
	}
	
	std::vector<reco::TransientTrack> transTrkVect = convToTransTrck(input);
	if(transTrkVect.size()!=3){
		if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: Bad transient track size = %d.\n", transTrkVect.size());
		return false;
	}
	TransientVertex secVtx;
	if(!checkSecVtx(transTrkVect, secVtx)){
		if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: Skip tauCand.\n");
		return false;
	}
	//	std::cout<<"sec vtx = "<<secVtx.position()<<std::endl;
	double massA1 = getInvariantMass(input, 0.140);//better construct from KinematicParticles?
	TLorentzVector lorentzA1 = getSumTLorentzVec(input, massA1);
	VertexRotation vtxC(lorentzA1, verbosity_);
	double thetaMax = vtxC.calcThetaMax();
	if(lorentzA1.M() > 2.0 || lorentzA1.M() < 3*.140){//soft upper value due to neutrino resolution
		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: Bad a1 mass = %f. Skip tauCand.\n", lorentzA1.M());
		return false;
	}
	
	double theta0;
	TVector3 tauFlghtDir;
	//bool worked = 
	vtxC.tryCorrection(primVtx, secVtx, theta0, tauFlghtDir);//can modify all 4 values; rotation is forced to reach thetaMax
	//	if(cnt == 26) vtxC.dumpEvt(primVtx, secVtx, thetaMax);
	//	if (!worked) return false;
	
	VirtualKinematicParticleFactory factory;
	for(unsigned int i = 0; i!=transTrkVect.size();i++){//apply vertex mod
		//		pions.push_back(kinFactory.particle(transTrkVect[i],piMass,piChi,piNdf,piMassSigma));
		RefCountedKinematicParticle tmp = kinFactory.particle(transTrkVect[i],piMass,piChi,piNdf,piMassSigma);
		//modify vertex
		AlgebraicVector7 newPar, oldPar = tmp->currentState().kinematicParameters().vector();
		newPar(0) = secVtx.position().x();
		newPar(1) = secVtx.position().y();
		newPar(2) = secVtx.position().z();
		newPar(3) = oldPar(3);
		newPar(4) = oldPar(4);
		newPar(5) = oldPar(5);
		newPar(6) = oldPar(6);
		KinematicParameters newParm(newPar);
		KinematicState newkineState(newParm, tmp->currentState().kinematicParametersError(), tmp->currentState().particleCharge(), transTrackBuilder_.field());
		float chi2 = tmp->chiSquared(), ndof = tmp->degreesOfFreedom();
		pions.push_back(factory.particle(newkineState, chi2, ndof, 0,0));
		if(verbosity_>=10){
			GlobalVector test = calcPVSVDir(primVtx, secVtx);
			TVector3 test2(test.x(),test.y(),test.z());
			double th = vtxC.unsignedAngle(test2, lorentzA1.Vect());
			printf("ThreeProngTauCreator::createStartScenario: delta GJ angle (%f).\n", th-thetaMax);
		}
	}
	
	if(theta0>TMath::Pi()/2){
		//		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: Unrealistic GJ angle (%f). But Event not skipped!?!\n", theta0);
		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: Unrealistic GJ angle (%f). tauCand skipped!\n", theta0);
		return false;
	}
	
	//??????always use thetaMax -> no amb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!and near to MC reality, but modification often out of vertex errors
	if ((thetaMax - theta0) < -pow(10.0,-10.0)*thetaMax && thetaMax>0) {
		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: thetaGJ = %f but replaced by thetaMax = %f.\n", theta0, thetaMax);
		//		if(!correctPrimVtx(primVtx, secVtx, tauFlghtDir, lorentzA1.Vect(), theta0, thetaMax)) return false;//could not correct PV within its errors
		theta0 = thetaMax;
	}
	
	//	else{//try to get rid of ambiguity
	//		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: thetaGJ = %f < thetaMax = %f. But try to use thetaMax.\n", theta0, thetaMax);
	//		if(correctPrimVtx(primVtx, secVtx, tauFlghtDir, lorentzA1.Vect(), theta0, thetaMax)) theta0 = thetaMax;
	//		//else keep theta = theta0
	//	}
	
	//	double thetaGen = genTau->momentum().Angle(genA1->momentum());
	//	std::cout<<" (thetaGen = "<<thetaGen<<")"<<std::endl;
	std::pair<double,double> tauSolutions = getTauMomentumMagnitudes(lorentzA1.M(),lorentzA1.P(),1.777,theta0);
	if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: tau solutions = %f, %f.\n", tauSolutions.first, tauSolutions.second);
	/*	if(tauSolutions.first==tauSolutions.second){
	 theta0 *= .5;
	 tauSolutions = getTauMomentumMagnitudes(lorentzA1.M(),lorentzA1.P(),1.777,theta0);
	 if(verbosity_>=2) printf("ThreeProngTauCreator::createStartScenario: thetaGJ is halved to %f, tau solutions are replaced by %f, %f.\n", theta0, tauSolutions.first, tauSolutions.second);
	 }*/
	bool ambiguity = false;
	if(fabs(tauSolutions.first-tauSolutions.second) < pow(10.0,-6)) ambiguity = true;
	//	std::pair<double,double> genTauSolutions = getTauMomentumMagnitudes(genA1->mass(),genA1->p(),1.777,thetaGen);
	//	std::cout<<"tau sols from gen: "<<genTauSolutions.first<<", "<<genTauSolutions.second<<std::endl;
	
	TLorentzVector tauGuess1;
	tauFlghtDir = tauFlghtDir.Unit();
	tauGuess1.SetXYZM(tauFlghtDir.X()*tauSolutions.first , tauFlghtDir.Y()*tauSolutions.first , tauFlghtDir.Z()*tauSolutions.first , 1.777);
	neutrinos.push_back(unknownNu(tauGuess1, lorentzA1, secVtx));
	if(ambiguity==true){
		TLorentzVector tauGuess2;
		tauGuess2.SetXYZM(tauFlghtDir.X()*tauSolutions.second, tauFlghtDir.Y()*tauSolutions.second, tauFlghtDir.Z()*tauSolutions.second, 1.777);
		neutrinos.push_back(unknownNu(tauGuess2, lorentzA1, secVtx));
	}
	if(neutrinos.size() != 1 && neutrinos.size() != 2){
		if(verbosity_>=1) printf("ThreeProngTauCreator::createStartScenario: Bad neutrino size = %d.\n", neutrinos.size());
		return false;
	}
	
	return true;
}
bool ThreeProngTauCreator::kinematicRefit(std::vector<RefCountedKinematicParticle> &unfitDaughters, const reco::Vertex &primVtx){
	if(unfitDaughters.size()!=4){
		printf("ThreeProngTauCreator::kinematicRefit:ERROR! Wrong size of daughters. Skip tauCand.\n");
		return false;
	}
	ParticleMass tauMass = 1.777;
	MultiTrackKinematicConstraint *tauMass_c = new  MultiTrackMassKinematicConstraint(tauMass, 4);
	
	std::vector<MultiTrackKinematicConstraint* > *constraintVector = new std::vector<MultiTrackKinematicConstraint* >;
	constraintVector->push_back(tauMass_c);
	GlobalPoint *primPoint = new GlobalPoint(primVtx.x(), primVtx.y(), primVtx.z());
	MultiTrackKinematicConstraint *pointing_c = new MultiTrackPointingKinematicConstraint(*primPoint);
	constraintVector->push_back(pointing_c);
	MultiTrackKinematicConstraint *combiC = new CombinedKinematicConstraint(*constraintVector);

	GlobalPoint vtxGuess = unfitDaughters[3]->currentState().globalPosition();//nu was created at common/corrected vertex of pions

	try{
		kinTree_ = kcvFitter_.fit(unfitDaughters, combiC, &vtxGuess);
	}catch(VertexException){//("KinematicStatePropagator without material::propagation failed!")
		std::cout<<"VertexException. Skip tauCand."<<std::endl;
		return false;
	}
	
	delete combiC;
	delete pointing_c;
	delete primPoint;
	delete constraintVector;
	delete tauMass_c;
	
	if(kinTree_->isValid()) return true;
	else{
		printf("ThreeProngTauCreator::kinematicRefit: ERROR! Tree is not valid. Skip tauCand.\n");
		return false;
	}
}

bool ThreeProngTauCreator::choose3bestTracks(std::vector<reco::TrackRef> &input, reco::Vertex & pVtx){
	int n = input.size(), m = 3;
	int noOfPerm = nOverM(n, m);
	if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks: find %d out of %d. need to compute %d permutations.\n", m, n, noOfPerm);
	sort(input.begin(), input.end(), cmpPt<reco::TrackRef>);
	std::vector<std::vector<reco::TrackRef> > combis = permuteCombinations(input);
	if(combis.size()==0){
		std::cout<<"ThreeProngTauCreator::choose3bestTracks:Bad combis. Skip it."<<std::endl;
		return false;
	}
	std::vector<std::pair<int,float> > chi2s;
	std::vector<std::pair<int,double> > movements;
	unsigned index=0;
	for (std::vector<std::vector<reco::TrackRef> >::iterator iter=combis.begin(); iter!=combis.end();) {
		if(verbosity_>=3) printf("pts: %f, %f, %f\n", iter->at(0)->pt(), iter->at(1)->pt(), iter->at(2)->pt());
		if(!sumCharge(*iter)){
			iter = combis.erase(iter);
			if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks: erased combi due to wrong charge sum. %d combis left.\n", combis.size());
			continue;
		}
		double massA1 = getInvariantMass(*iter, 0.140);
		if(massA1 > 2.0 || massA1 < 3*0.140){//soft upper value
			iter = combis.erase(iter);
			if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks: erased combi due to wrong mass. %d combis left.\n", combis.size());
			continue;
		}
		TransientVertex tmpVtx;
		std::vector<reco::TransientTrack> trks = convToTransTrck(*iter);
		if(!checkSecVtx(trks, tmpVtx)){
			iter = combis.erase(iter);
			if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks: erased combi due to bad vertex. %d combis left.\n", combis.size());
			continue;
		}
		//			printf("vtx = (%8.6f, %8.6f, %8.6f), chi2ndof = %f\n", tmpVtx.position().x(), tmpVtx.position().y(), tmpVtx.position().z(), tmpVtx.normalisedChiSquared());
		
		TLorentzVector lorentzA1 = getSumTLorentzVec(*iter, massA1);
		VertexRotation vtxC(lorentzA1, verbosity_);
		double theta0;
		TVector3 tauFlghtDir;
		reco::Vertex pvTemp = pVtx;//do not modify original pv here
		vtxC.tryCorrection(pvTemp, tmpVtx, theta0, tauFlghtDir, false);//do not force rotation, only rotate within errors
		if(vtxC.isValid()) movements.push_back(std::make_pair(index,vtxC.movement()));
		else{
			iter = combis.erase(iter);
			if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks: erased combi due to bad vertex correction. %d combis left.\n", combis.size());
			continue;
		}
		
		chi2s.push_back(std::make_pair(index,tmpVtx.normalisedChiSquared()));
		
		++index;
		++iter;//only moved if nothing was deleted
	}
	if (combis.size()<1){
		if (combis.size()>=1) printf("ThreeProngTauCreator::choose3bestTracks:No combi survived.\n");
		return false;
	}
	
	if (combis.size()>1){//chose the pair with smallest vertex rotation needed!!!
		//		sort(chi2s.begin(), chi2s.end(), pairSecond<int, float>);
		//		if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks:Too much combis (%d) left. Take best common vtx (chi2/ndf = %f), second best is (chi2/ndf = %f).\n", combis.size(), chi2s.at(0).second, chi2s.at(1).second);
		//		unsigned int i = chi2s.front().first;
		//		input.assign( combis.at(i).begin(), combis.at(i).end() );
		
		sort(movements.begin(), movements.end(), pairSecond<int, double>);
		if(verbosity_>=2) printf("ThreeProngTauCreator::choose3bestTracks:Too much combis (%d) left. Take one with smallest vtx correction movement (%f [sigma]), second best is (%f [sigma]).\n", combis.size(), movements.front().second, movements.at(1).second);
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
	//TransientTrackBuilder liefert reco::TransientTrack*
	std::vector<reco::TransientTrack> transTrkVct;
	for (std::vector<reco::TrackRef>::iterator iter=input.begin(); iter!=input.end(); ++iter) {
		transTrkVct.push_back( transTrackBuilder_.build( *iter ) );
	}
	return transTrkVct;
}
bool ThreeProngTauCreator::checkSecVtx(std::vector<reco::TransientTrack> &trkVct, TransientVertex & transVtx){
	if(trkVct.size()<2){
		printf("Can't check SecVertex: Only %i Tracks.", trkVct.size());
		return false;
	}else{
		bool useAdaptive = false;
		if(useAdaptive){
			AdaptiveVertexFitter avf;//more robust
			avf.setWeightThreshold(0.1);//weight per track. lasse (fast) alle fits zu, sonst exception
			transVtx = avf.vertex(trkVct);//AdaptiveVertexFitter
		}else{
			KalmanVertexFitter kvf(true);
			transVtx = kvf.vertex(trkVct);//KalmanVertexFitter
		}
		//("TwoTrackMinimumDistance ... comparing track with itself!")
		if(!transVtx.isValid()) printf("ThreeProngTauCreator::checkSecVtx: Secondary vertex not valid.\n");
		else if(verbosity_>=2) printf("ThreeProngTauCreator::checkSecVtx: vtx(%f, %f, %f)pm(%f, %f, %f).\n", transVtx.position().x(), transVtx.position().y(), transVtx.position().z(), transVtx.positionError().cxx(), transVtx.positionError().cyy(), transVtx.positionError().czz());
		if(!useAdaptive){
			if(!transVtx.hasRefittedTracks()){
				if(verbosity_>=2) printf("ThreeProngTauCreator::checkSecVtx: Secondary has 0 refitted tracks.\n");
				return false;
			}else if(transVtx.refittedTracks().size()!=trkVct.size()){
				if(verbosity_>=2) printf("ThreeProngTauCreator::checkSecVtx: Secondary has only %i refitted of %i initial tracks.\n", transVtx.refittedTracks().size(), trkVct.size());
				return false;
			}
		}
		
		return transVtx.isValid();
	}
}
GlobalVector ThreeProngTauCreator::calcPVSVDir(reco::Vertex &primVtx, TransientVertex & transVtx){//calculate direction from primVtx to secVtx
	GlobalPoint *secVtx = new GlobalPoint(transVtx.position());
	GlobalVector primVertex(primVtx.x(),primVtx.y(),primVtx.z());
	
	GlobalVector dirPSVS(secVtx->x(),secVtx->y(),secVtx->z());
	dirPSVS-=primVertex;
	//		std::cout<<"unitVct "<<dirPSVS.unit()<<std::endl;
	
	//		return dirPSVS.unit();
	return dirPSVS;//do not use unit as one will try to rotate pv around sv
}
std::pair<double,double> ThreeProngTauCreator::getTauMomentumMagnitudes(double ma1,double pa1,double M,double theta){
	//	std::cout<<"\n(ma1,pa1,theta) "<<ma1<<","<<pa1<<","<<theta<<std::endl;
	double root = sqrt((pow(ma1,2) + pow(pa1,2))*(pow(pow(ma1,2) - pow(M,2),2) -4*pow(M,2)*pow(pa1,2)*pow(sin(theta),2))),
	numerator = (pow(ma1,2) + pow(M,2))*pa1*cos(theta),
	denominator = (2*pow(ma1,2) + 2*pow(pa1,2)*pow(sin(theta),2));
	double tauMomentumMagnitude1 = (numerator - root) / denominator,
	tauMomentumMagnitude2 = (numerator + root)	/ denominator;
	
	if(!(tauMomentumMagnitude1>0)) tauMomentumMagnitude1 = numerator / denominator;//fange 'nan' ab und ersetze durch Quasi-Mittelwert der beiden solutions
	if(!(tauMomentumMagnitude2>0)) tauMomentumMagnitude2 = tauMomentumMagnitude1;
	
	if(!(tauMomentumMagnitude1>0)) tauMomentumMagnitude1 = 0;//fange 'nan' ab
	if(!(tauMomentumMagnitude2>0)) tauMomentumMagnitude2 = 0;
	
	return std::make_pair(tauMomentumMagnitude1,tauMomentumMagnitude2);
}
RefCountedKinematicParticle ThreeProngTauCreator::virtualKinematicParticle(GlobalPoint vtxGuess, GlobalVector impulsGuess){
	VirtualKinematicParticleFactory factory;
	//(x,y,z,p_x,p_y,p_z,m)
	const KinematicParameters parameters(AlgebraicVector7(vtxGuess.x(),vtxGuess.y(),vtxGuess.z(),impulsGuess.x(),impulsGuess.y(),impulsGuess.z(),pow(10.,-10.)));//start-pT aus MET nehmen?
	//SymMatrix constructed out of lower/upper (default lower) block data with dimension N*(N+1)/2
	ROOT::Math::SVector<double,28> svector28;//(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,100);//dxm und co sind 0, da unkorreliert
	for(unsigned int i=1; i!=22; i++) svector28(i-1) = pow(10.,-12.);//0.;//korrelationen
	for(unsigned int i=22; i!=28; i++) svector28(i-1) = pow(10.,-12.);//0.;//korrelationen zw m und p/vertex
	for(unsigned int n=1; n!=7; n++) svector28(n*(n+1)/2 - 1) = pow(10.,2.);//Diagonalelemente, HugeError-Methode
	svector28(27) = pow(10.,-12.);//massen-fehler
	
	// ROOT::Math::SVectorr<double, 6> v(1,2,3,4,5,6);
	//	SMatrixSym3 s1(v);  // lower block (default)
	// this will produce the symmetric  matrix
	//    (  1    2    4 
	//       2    3    5
	//       4    5    6  )
	//	SMatrixSym3 s2(v,false);  // upper block 
	// this will produce the symmetric  matrix
	//    (  1    2    3
	//       2    4    5
	//       3    5    6  )
	//	AlgebraicSymMatrix77 matrix(svector28);
	ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > matrix(svector28);
	//	matrix.Print(std::cout);printf("\n");
	
	const KinematicParametersError parametersError(matrix);
	const TrackCharge charge = 0;
	KinematicState kineState(parameters, parametersError, charge, transTrackBuilder_.field());
	float chiSquared=1000, degreesOfFr=7;
	
	return factory.particle(kineState, chiSquared, degreesOfFr, 0,0);
}
RefCountedKinematicParticle ThreeProngTauCreator::unknownNu(TLorentzVector &tauGuess, TLorentzVector &a1, TransientVertex & secVtx){
	GlobalVector nuImpGuess(1,1,1);
	GlobalVector tauImpGuess(tauGuess.X(), tauGuess.Y(), tauGuess.Z());
	
	TLorentzVector nuGuess = tauGuess-a1;
	if(verbosity_>=1) std::cout<<"ThreeProngTauCreator::unknownNu: nuGuess (px,py,pz,m) "<<nuGuess.Px()<<","<<nuGuess.Py()<<","<<nuGuess.Pz()<<","<<nuGuess.M()<<std::endl;
	//	TLorentzVector nuGuess1(0,0,0,0);
	//	nuGuess1.SetVect(tauGuess1.Vect() - recoA1->momentum());
	//	std::cout<<"nuGuess1.P()="<<nuGuess1.P()<<", (px="<<nuGuess1.X()<<", py="<<nuGuess1.Y()<<", pz="<<nuGuess1.Z()<<")"<<std::endl;
	if(tauGuess.P()==0)  nuGuess.SetXYZM(1,1,1,0);
	nuImpGuess = GlobalVector(nuGuess.X(),nuGuess.Y(),nuGuess.Z());
	return virtualKinematicParticle(secVtx.position(), nuImpGuess);
}
template <typename T> std::vector<std::vector<T> > ThreeProngTauCreator::permuteCombinations(const std::vector<T> &vect){
	std::vector<std::vector<T> > combis;
	int m = 3, n = vect.size();
	int unsigned combinations = nOverM(n, m);
	//		printf("calc %d combis\n", combinations);
//	int cnt = 0;
	typename std::vector<T>::const_iterator iter1, iter2, iter3;
	for (iter1 = vect.begin(); iter1!=vect.end(); ++iter1) {
		//			printf("iter1 %d\n", *iter1);
		iter2 = iter1;
		++iter2;
		for (; iter2!=vect.end(); ++iter2) {
			//				printf("iter2 %d\n", *iter2);
			iter3 = iter2;
			++iter3;
			for (; iter3!=vect.end(); ++iter3) {
				//					printf("iter3 %d\n", *iter3);
//				cnt++;
				std::vector<T> newCombi;
				newCombi.push_back(*iter1);
				newCombi.push_back(*iter2);
				newCombi.push_back(*iter3);
				combis.push_back(newCombi);
				//					printf("combi %2d: %d, %d, %d\n", cnt, *iter1, *iter2, *iter3);
			}			
		}
	}
	if(combis.size() != combinations) combis.clear();
	return combis;
}
