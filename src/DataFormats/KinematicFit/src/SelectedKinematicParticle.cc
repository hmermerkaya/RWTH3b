#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"

SelectedKinematicParticle::SelectedKinematicParticle() {
  status_ = -1;
  name_ = "";
  charge_ = 0;
  ambiguity_ = -1;
  
  kinparm_.ResizeTo(7);
  input_kinparm_.ResizeTo(7);
  kinmatrix_.ResizeTo(TMatrixDSym(7));
  input_kinmatrix_.ResizeTo(TMatrixDSym(7));
  
  CandRef_ = reco::RecoChargedCandidateRef();
}
SelectedKinematicParticle::SelectedKinematicParticle(const RefCountedKinematicParticle & kinparticle, const int status, const std::string & name, const int ambiguity, const reco::RecoChargedCandidateRef & CandRef){
  status_ = status;
  name_ = name;
  charge_ = (kinparticle->currentState()).particleCharge();
  ambiguity_ = ambiguity;
  
  kinparm_.ResizeTo(7);
  kinparm_ = convertVector((kinparticle->currentState()).kinematicParameters().vector());
  input_kinparm_.ResizeTo(7);
  input_kinparm_ = convertVector((kinparticle->initialState()).kinematicParameters().vector());
  
  kinmatrix_.ResizeTo(TMatrixDSym(7));
  kinmatrix_ = convertMatrix((kinparticle->currentState()).kinematicParametersError().matrix());
  input_kinmatrix_.ResizeTo(TMatrixDSym(7));
  input_kinmatrix_ = convertMatrix((kinparticle->initialState()).kinematicParametersError().matrix());
  
  CandRef_ = CandRef;
}

const int SelectedKinematicParticle::status() const {
  return status_;
}
const std::string & SelectedKinematicParticle::name() const {
  return name_;
}
const int SelectedKinematicParticle::charge() const {
  return charge_;
}
const int SelectedKinematicParticle::ambiguity() const {
  return ambiguity_;
}

const TVectorT<double> & SelectedKinematicParticle::parameters() const {
  return kinparm_;
}
const TVectorT<double> & SelectedKinematicParticle::input_parameters() const {
  return input_kinparm_;
}

const TMatrixDSym & SelectedKinematicParticle::matrix() const {
  return kinmatrix_;
}
const TMatrixDSym & SelectedKinematicParticle::input_matrix() const {
  return input_kinmatrix_;
}

const reco::RecoChargedCandidateRef & SelectedKinematicParticle::candRef() const {
  return CandRef_;
}
void SelectedKinematicParticle::setCandRef(const reco::RecoChargedCandidateRef & parm){
  CandRef_ = parm;
}

const TVectorT<double> SelectedKinematicParticle::convertVector( const AlgebraicVector7 & vector ) {
  TVectorT<double> tmpVector(7);
  for (int i = 0; i < 7; i++) {
    tmpVector[i] = double(vector.At(i));
  }
  return tmpVector;
}
const TMatrixDSym SelectedKinematicParticle::convertMatrix( const AlgebraicSymMatrix77 & matrix ) {
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

const TLorentzVector SelectedKinematicParticle::p4() const {
  TLorentzVector p4tmp;
  p4tmp.SetVectM(TVector3(kinparm_[3], kinparm_[4], kinparm_[5]), kinparm_[6]);
  return p4tmp;
}
const TVector3 SelectedKinematicParticle::vertex() const {
  TVector3 vtx(kinparm_[0], kinparm_[1], kinparm_[2]);
  return vtx;
}

void SelectedKinematicParticle::setInitialState(const TLorentzVector & momentum, const reco::Vertex & primVtx) {
  input_kinparm_.ResizeTo(7);
  input_kinparm_[0] = primVtx.x();
  input_kinparm_[1] = primVtx.y();
  input_kinparm_[2] = primVtx.z();
  input_kinparm_[3] = momentum.Px();
  input_kinparm_[4] = momentum.Py();
  input_kinparm_[5] = momentum.Pz();
  input_kinparm_[6] = momentum.M();
  input_kinmatrix_.ResizeTo(TMatrixDSym(7));
  for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	input_kinmatrix_[i][j] = primVtx.covariance(i,j);//how to treat correlation between vertex and momentum block?
}

