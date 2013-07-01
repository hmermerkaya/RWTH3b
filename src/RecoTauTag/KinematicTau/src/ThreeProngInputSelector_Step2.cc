#include "RecoTauTag/KinematicTau/interface/ThreeProngInputSelector_Step2.h"
#include "RecoTauTag/KinematicTau/interface/SecondaryVertexHelper.h"

ThreeProngInputSelector_Step2::ThreeProngInputSelector_Step2(const edm::ParameterSet & iConfig):
  primVtxTag_(iConfig.getParameter<edm::InputTag>("primVtx")),
  KinematicTauCandTag_(iConfig.getParameter<edm::InputTag>("KinematicTauCandTag")),
  VertexTags_(iConfig.getUntrackedParameter< std::vector<std::string> >("VertexTags")),
  TauVtxList_(iConfig.getUntrackedParameter< std::vector<std::string> >("NonTauTracks")),
  minTau_(iConfig.getUntrackedParameter<unsigned int>("minTau", 1)), //filter returns true if more/equal than minTau_ taus were selected
  etacut_(iConfig.getUntrackedParameter<double>("etacut",2.1)),
  sigcut_(iConfig.getUntrackedParameter<double>("sigcut",3.0))
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
  
  std::auto_ptr<std::vector<SelectedKinematicDecay> > PreKinematicDecaysStep2 =  std::auto_ptr<std::vector<SelectedKinematicDecay> >(new std::vector<SelectedKinematicDecay>);
  std::vector<SelectedKinematicDecay> &PreKinematicDecaysStep2_=*PreKinematicDecaysStep2;

  select(PreKinematicDecaysStep2_,iSetup);
  
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

bool ThreeProngInputSelector_Step2::select(std::vector<SelectedKinematicDecay> &PreKinematicDecaysStep2_,const edm::EventSetup& iSetup) {
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder_;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder_);

  //load primary vertex collection created w/o tau tracks (needs to be moved later)
  edm::Handle<reco::VertexCollection > primaryVertexCollection;
  iEvent_->getByLabel(primVtxTag_, primaryVertexCollection);
  if (!primaryVertexCollection.isValid()) return false;
  reco::VertexCollection primaryVertices = *primaryVertexCollection;
  if(primaryVertices.size()>0){
    
    edm::Handle<std::vector<std::vector<SelectedKinematicDecay> > > KinematicTauCandidate;
    iEvent_->getByLabel(KinematicTauCandTag_,KinematicTauCandidate);
    
    for(unsigned int i=0; i<KinematicTauCandidate->size();i++){
      SelectedKinematicDecay bestTau;
      bool hasTau=false;
      for(unsigned int j=0; j<KinematicTauCandidate->at(i).size();j++){
	SelectedKinematicDecay KTau=KinematicTauCandidate->at(i).at(j);
	if(j==0){
	  const reco::PFTauRef pfiter=KTau.PFTauRef();
	}
	SecondaryVertexHelper SVH(transTrackBuilder_,KTau);
	if(SVH.hasSecondaryVertex()){
	  TString vertexName=KTau.PrimaryVertexReFitCollectionTag();
	  TString VTag;
	  for(unsigned int v=0;v<VertexTags_.size() && VertexTags_.size()==TauVtxList_.size();v++){
	    if(vertexName==TauVtxList_.at(v)) VTag=VertexTags_.at(v);
	  }
	  edm::Handle<reco::VertexCollection > CurrentTauPrimaryVtx;
	  iEvent_->getByLabel(edm::InputTag(VTag.Data()),CurrentTauPrimaryVtx);
	  if(!CurrentTauPrimaryVtx.isValid()) continue;
	  if(CurrentTauPrimaryVtx->size()==0) continue;
	  reco::Vertex primaryVertexReFit=CurrentTauPrimaryVtx->front();
	  reco::Vertex primaryVertexReFitAndRotated=primaryVertexReFit;
	  TVector3 tauFlghtDirNoCorr;
	  TVector3 tauFlghtDir;
	  TLorentzVector a1_p4=SVH.Initial_a1_p4();
	  double initThetaGJ,ThetaMax;
	  TransientVertex SecondaryVertex=SVH.InitialSecondaryVertex();
	  std::vector<reco::TransientTrack> RefittedTracks=SVH.InitialRefittedTracks();
	  double s = VertexRotationAndSignificance(SecondaryVertex,RefittedTracks,tauFlghtDirNoCorr,
						   primaryVertexReFitAndRotated,a1_p4,tauFlghtDir,initThetaGJ,ThetaMax);
	  if(/*(s<sigcut_  || fabs(initThetaGJ)<fabs(ThetaMax)) &&*/ s>=0 ){// prevent nan
	    if(fabs(tauFlghtDirNoCorr.Eta())<etacut_ || fabs(tauFlghtDir.Eta())<etacut_ ){
	      KTau.SetInitialVertexProperties(primaryVertexReFit,primaryVertexReFitAndRotated,
					      SVH.InitialRefittedTracks(),SVH.InitialSecondaryVertex());
	      KTau.SetInitialKinematics(tauFlghtDirNoCorr.Unit(),SVH.Initial_pions(),a1_p4,tauFlghtDir.Unit(),initThetaGJ,ThetaMax);
	      
	      bestTau=KTau;
	      hasTau=true;
	    }
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

double ThreeProngInputSelector_Step2::VertexRotationAndSignificance(TransientVertex &tmpVtx, std::vector<reco::TransientTrack> trks, 
								    TVector3 &tauFlghtDirNoCorr,
								    reco::Vertex &pVtx, TLorentzVector &lorentzA1, 
								    TVector3 &tauFlghtDir,double &theta0, double &thetaMax){
  TVector3 pv(pVtx.position().x(), pVtx.position().y(), pVtx.position().z());
  TVector3 sv(tmpVtx.position().x(), tmpVtx.position().y(), tmpVtx.position().z());
  tauFlghtDirNoCorr = sv - pv;

  VertexRotation vtxC(lorentzA1);
  thetaMax=fabs(vtxC.calcThetaMax());
  return vtxC.rotatePV(pVtx,tmpVtx,theta0, tauFlghtDir);
}


//define this as a plug-in
DEFINE_FWK_MODULE(ThreeProngInputSelector_Step2);
