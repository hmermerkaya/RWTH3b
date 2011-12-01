#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"

SelectedKinematicParticle::SelectedKinematicParticle() {
    status__ = -1;
    matched__ = -1;

    kinparm__.ResizeTo(7);
    input_kinparm__.ResizeTo(7);
    kinmatrix__.ResizeTo(TMatrixDSym(7));
    input_kinmatrix__.ResizeTo(TMatrixDSym(7));
}
SelectedKinematicParticle::SelectedKinematicParticle(const TVectorT<double> kinparm_, const TMatrixDSym kinmatrix_, const int charge_, const std::string name_, const float chi2_, const float ndf_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_ = 0, const int status_ = 0) {
    status__ = status_;
    matched__ = -1;
    iterations__ = iterations_;
    maxiterations__ = maxiterations_;
    ambiguity__ = ambiguity_;
    
    charge__ = charge_;
    chi2__ = chi2_; 
	ndf__ = ndf_;
    csum__ = csum_;
    mincsum__ = mincsum_;
    
    name__ = name_;

    kinparm__.ResizeTo(7);
	kinparm__ = kinparm_;
    input_kinparm__.ResizeTo(7);
    
    kinmatrix__.ResizeTo(TMatrixDSym(7));
	kinmatrix__ = kinmatrix_;
    input_kinmatrix__.ResizeTo(TMatrixDSym(7));
    
    CandRef__ = CandRef_;
}
SelectedKinematicParticle::SelectedKinematicParticle(const TVectorT<double> kinparm_, const TMatrixDSym kinmatrix_, const TVectorT<double> input_kinparm_, const TMatrixDSym input_kinmatrix_, const int charge_, const std::string name_, const float chi2_, const float ndf_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_ = 0, const int status_ = 0){
    status__ = status_;
    matched__ = -1;
    iterations__ = iterations_;
    maxiterations__ = maxiterations_;
    ambiguity__ = ambiguity_;
    
    charge__ = charge_;
    chi2__ = chi2_; 
	ndf__ = ndf_;
    csum__ = csum_;
    mincsum__ = mincsum_;
    
    name__ = name_;
    
    kinparm__.ResizeTo(7);
	kinparm__ = kinparm_;
    input_kinparm__.ResizeTo(7);
    input_kinparm__ = input_kinparm_;
    
    kinmatrix__.ResizeTo(TMatrixDSym(7));
	kinmatrix__ = kinmatrix_;
    input_kinmatrix__.ResizeTo(TMatrixDSym(7));
    input_kinmatrix__ = input_kinmatrix_;
    
    CandRef__ = CandRef_;    
}
SelectedKinematicParticle::SelectedKinematicParticle(const RefCountedKinematicParticle kinparticle_, const std::string name_, const int iterations_, const int maxiterations_, const float csum_, const float mincsum_, const reco::RecoChargedCandidateRef CandRef_, const int ambiguity_, const int status_){
    status__ = status_;
    matched__ = -1;
    iterations__ = iterations_;
    maxiterations__ = maxiterations_;
    ambiguity__ = ambiguity_;
    
    charge__ = (kinparticle_->currentState()).particleCharge();
    chi2__ = kinparticle_->chiSquared(); 
	ndf__ = kinparticle_->degreesOfFreedom();
    csum__ = csum_;
    mincsum__ = mincsum_;
    
    name__ = name_;
    
    kinparm__.ResizeTo(7);
	kinparm__ = convertVector((kinparticle_->currentState()).kinematicParameters().vector());
    input_kinparm__.ResizeTo(7);
    input_kinparm__ = convertVector((kinparticle_->initialState()).kinematicParameters().vector());
    
    kinmatrix__.ResizeTo(TMatrixDSym(7));
	kinmatrix__ = convertMatrix((kinparticle_->currentState()).kinematicParametersError().matrix());
    input_kinmatrix__.ResizeTo(TMatrixDSym(7));
    input_kinmatrix__ = convertMatrix((kinparticle_->initialState()).kinematicParametersError().matrix());
    
    CandRef__ = CandRef_;    
    
}
int SelectedKinematicParticle::status() const {
    return status__;
}
int SelectedKinematicParticle::matched() const {
    return matched__;
}
int SelectedKinematicParticle::iterations() const {
    return iterations__;
}
int SelectedKinematicParticle::maxiterations() const {
    return maxiterations__;
}
int SelectedKinematicParticle::ambiguity() const {
    return ambiguity__;
}
int SelectedKinematicParticle::charge() const {
    return charge__;
}
float SelectedKinematicParticle::chi2() const {
    return chi2__;
}
float SelectedKinematicParticle::ndf() const {
    return ndf__;
}
float SelectedKinematicParticle::csum() const {
    return csum__;
}
float SelectedKinematicParticle::mincsum() const {
    return mincsum__;
}
std::string SelectedKinematicParticle::name() const {
    return name__;
}
TVectorT<double> SelectedKinematicParticle::parameters() const {
    return kinparm__;
}
TVectorT<double> SelectedKinematicParticle::input_parameters() const {
    return input_kinparm__;
}
TMatrixDSym SelectedKinematicParticle::matrix() const {
    return kinmatrix__;
}
TMatrixDSym SelectedKinematicParticle::input_matrix() const {
    return input_kinmatrix__;
}
reco::RecoChargedCandidateRef SelectedKinematicParticle::candRef() const {
    return CandRef__;
}
void SelectedKinematicParticle::setCandRef(const reco::RecoChargedCandidateRef parm){
	CandRef__ = parm;
}

TVectorT<double> SelectedKinematicParticle::convertVector( const AlgebraicVector7 vector ) {
	TVectorT<double> tmpVector(7);
	for (int i = 0; i < 7; i++) {
		tmpVector[i] = double(vector.At(i));
	}
	return tmpVector;
}
TMatrixDSym SelectedKinematicParticle::convertMatrix( const AlgebraicSymMatrix77 matrix ) {
	TMatrixDSym tmpMatrix(7);
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			tmpMatrix[i][j] = matrix.At(i,j);
            if (i == j && tmpMatrix[i][j] < -0.000001) {
                //std::cerr << "WARNING: matrix_el[" << i << "][" << j << "]: " << tmpMatrix[i][j] << " (before: " << matrix.At(i,j) << ")" << std::endl;
            }
		}
	}
	return tmpMatrix;
}

TLorentzVector SelectedKinematicParticle::p4() const {
	TLorentzVector p4tmp;
	p4tmp.SetVectM(TVector3(kinparm__[3], kinparm__[4], kinparm__[5]), kinparm__[6]);
	return p4tmp;
}
TVector3 SelectedKinematicParticle::vertex() const {
	TVector3 vtx(kinparm__[0], kinparm__[1], kinparm__[2]);
	return vtx;
}

void SelectedKinematicParticle::setMatched(const int parm) {
    matched__ = parm;
}

void SelectedKinematicParticle::setInitialState(const TLorentzVector & momentum, const reco::Vertex & primVtx) {
    input_kinparm__.ResizeTo(7);
	input_kinparm__[0] = primVtx.x();
	input_kinparm__[1] = primVtx.y();
	input_kinparm__[2] = primVtx.z();
	input_kinparm__[3] = momentum.Px();
	input_kinparm__[4] = momentum.Py();
	input_kinparm__[5] = momentum.Pz();
	input_kinparm__[6] = momentum.M();

	input_kinmatrix__.ResizeTo(TMatrixDSym(7));
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	input_kinmatrix__[i][j] = primVtx.covariance(i,j);//how to treat correlation between vertex and momentum block?
}
	
