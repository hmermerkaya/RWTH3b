#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step2.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"

ThreeProngInputSelector_Step2::ThreeProngInputSelector_Step2(const edm::ParameterSet & iConfig):
  KinematicTauTools(),
  primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag")),
  minTau_(iConfig.getUntrackedParameter<unsigned int>("minTau", 1)) //filter returns true if more/equal than minTau_ taus were selected
{
    iConfig_ = iConfig;
    produces<std::vector<SelectedKinematicDecay> >("PreKinematicDecaysStep2");
}


ThreeProngInputSelector_Step2::~ThreeProngInputSelector_Step2() {
}

// ------------ method called on each new Event  ------------
void ThreeProngInputSelector_Step2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  cnt_++;
  iEvent_ = &iEvent;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);
  Set_TransientTrackBuilder(transTrackBuilder_);
  
  std::auto_ptr<std::vector<SelectedKinematicDecay> > PreKinematicDecaysStep2 =  std::auto_ptr<std::vector<SelectedKinematicDecay> >(new std::vector<SelectedKinematicDecay>);
  std::vector<SelectedKinematicDecay> &PreKinematicDecaysStep2_=*PreKinematicDecaysStep2;

  select(PreKinematicDecaysStep2_);
  
  iEvent_->put(PreKinematicDecaysStep2,"PreKinematicDecaysStep2");
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeProngInputSelector_Step2::beginJob() {
  cnt_ = 0;
  cntFound_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeProngInputSelector_Step2::endJob() {
  float ratio = 0.0;
  if (cnt_ != 0) ratio = (float)cntFound_/cnt_;
  edm::LogVerbatim("ThreeProngInputSelector_Step2") << "--> [ThreeProngInputSelector_Step2] found at least " << minTau_ << " 3-prongs per event. Efficiency: " << cntFound_ << "/" << cnt_ << " = " << std::setprecision(4) << ratio*100.0 << "%";
}

bool ThreeProngInputSelector_Step2::select(std::vector<SelectedKinematicDecay> &PreKinematicDecaysStep2_) {
  //load primary vertex collection created w/o tau tracks (needs to be moved later)
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtxTag_, primaryVertexCollection);
  if (!primaryVertexCollection.isValid()) return false;
  reco::VertexCollection primaryVertices = *primaryVertexCollection;
  if(primaryVertices.size()>0){
    
    edm::Handle<std::vector<std::vector<SelectedKinematicDecay> > > KinematicTauCandidate;
    iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidate);
    
    for(unsigned int i=0; i<KinematicTauCandidate->size();i++){
      double sig=1e13;// use high value to prevent missing tau
      SelectedKinematicDecay bestTau;
      bool hasTau=false;
      edm::LogInfo("ThreeProngInputSelector_Step2")<<"InputTrackSelector::select: i = " << i;
      for(unsigned int j=0; j<KinematicTauCandidate->at(i).size();j++){
	SelectedKinematicDecay KTau=KinematicTauCandidate->at(i).at(j);
	if(j==0){
	  const reco::PFTauRef pfiter=KTau.PFTauRef();
	  edm::LogInfo("ThreeProngInputSelector_Step2")<<"ThreeProngInputSelector_Step2::select:  PFTau "<< i <<" (Px,Py,Pz,E)=("<<(pfiter)->p4().Px() << ","<<(pfiter)->p4().Px() << ","<<(pfiter)->p4().Pz() << ","<<(pfiter)->p4().E() << ")";
	}
	SecondaryVertexHelper SVH(transTrackBuilder_,KTau);
	if(SVH.hasSecondaryVertex()){
	  reco::Vertex primaryVertexReFit=primaryVertices.front();
	  reco::Vertex primaryVertexReFitAndRotated=primaryVertexReFit;
	  TVector3 tauFlghtDir;
	  TLorentzVector a1_p4;
	  double initThetaGJ,ThetaMax;
	  TransientVertex SecondaryVertex=SVH.InitalSecondaryVertex();
	  std::vector<reco::TransientTrack> RefittedTracks=SVH.InitalRefittedTracks();
	  double s = VertexRotationAndSignificance(KTau.InitalTrackTriplet(),SecondaryVertex,RefittedTracks,primaryVertexReFitAndRotated,a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	  edm::LogInfo("ThreeProngInputSelector_Step2")<<"ThreeProngInputSelector_Step2::select: significance " << s 
						       << " PVertexFit and Rotate (" <<  primaryVertexReFitAndRotated.position().x() 
						       << "," <<  primaryVertexReFitAndRotated.position().y() 
						       << "," <<  primaryVertexReFitAndRotated.position().z() << ")" 
						       << " SV (" <<  SVH.InitalSecondaryVertex().position().x()
						       << "," <<  SVH.InitalSecondaryVertex().position().y()
						       << "," <<  SVH.InitalSecondaryVertex().position().z() << ")";
	  
	  if(s<sig && s>=0){// prevent nan
	    KTau.SetInitalVertexProperties(primaryVertexReFit,primaryVertexReFitAndRotated,SVH.InitalRefittedTracks(),SVH.InitalSecondaryVertex());
	    KTau.SetInitalKinematics(tauFlghtDir,SVH.Inital_pions(),a1_p4,initThetaGJ,ThetaMax);
	    bestTau=KTau;
	    hasTau=true;
	  }
	}
      }
      if(hasTau) PreKinematicDecaysStep2_.push_back(bestTau);
    }
    if (PreKinematicDecaysStep2_.size() >= minTau_) {
      LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: " << PreKinematicDecaysStep2_.size() << " tau candidate(s) reconstructed.";
      cntFound_++;
      return true;
    } 
    else {
      LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Warning! Only " << PreKinematicDecaysStep2_.size() << " tau candidate(s) reconstructed. Skip Event.";
      return false;
    }
  }
  LogTrace("ThreeProngInputSelector_Step2") << "ThreeProngInputSelector_Step2::select: Unable to calculate a new primary vertex. No tau candidates reconstructed.";
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step2);
