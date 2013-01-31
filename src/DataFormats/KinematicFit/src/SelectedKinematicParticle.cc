#include "DataFormats/KinematicFit/interface/SelectedKinematicParticle.h"

SelectedKinematicParticle::SelectedKinematicParticle() {
  status_ = -1;
  pdgid_=0;
  charge_ = 0;
  ambiguity_ = -1;
  
  kinparm_.ResizeTo(7);
  input_kinparm_.ResizeTo(7);
  kinmatrix_.ResizeTo(TMatrixDSym(7));
  input_kinmatrix_.ResizeTo(TMatrixDSym(7));
  
  CandRef_ = reco::RecoChargedCandidateRef();
}

SelectedKinematicParticle::SelectedKinematicParticle(LorentzVectorParticle p, const int status, const int ambiguity, const reco::RecoChargedCandidateRef & CandRef){
  status_ = status;
  pdgid_=p.PDGID();
  charge_ = p.Charge();
  ambiguity_ = ambiguity;
  
  kinparm_.ResizeTo(LorentzVectorParticle::NLorentzandVertexPar);
  for(int i=0;i<LorentzVectorParticle::NLorentzandVertexPar;i++) kinparm_(i) = p.Parameter(i);
  input_kinparm_.ResizeTo(LorentzVectorParticle::NLorentzandVertexPar);
  input_kinparm_ = kinparm_;
  
  kinmatrix_.ResizeTo(LorentzVectorParticle::NLorentzandVertexPar,LorentzVectorParticle::NLorentzandVertexPar);
  for(int i=0;i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    for(int j=0;j<LorentzVectorParticle::NLorentzandVertexPar;j++){kinmatrix_ = p.Covariance(i,j);}
  }
  input_kinmatrix_.ResizeTo(LorentzVectorParticle::NLorentzandVertexPar,LorentzVectorParticle::NLorentzandVertexPar);
  input_kinmatrix_ = kinmatrix_;
  
  CandRef_ = CandRef;
}

const int SelectedKinematicParticle::status() const {
  return status_;
}
const int  SelectedKinematicParticle::pdgid() const {
  return pdgid_;
}
const int SelectedKinematicParticle::charge() const {
  return charge_;
}
const unsigned int SelectedKinematicParticle::ambiguity() const {
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

