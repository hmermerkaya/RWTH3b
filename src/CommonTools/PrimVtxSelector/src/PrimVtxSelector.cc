#include "../interface/PrimVtxSelector.h"

PrimVtxSelector::PrimVtxSelector(const edm::ParameterSet& iConfig):
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
minTracks_(iConfig.getUntrackedParameter("minTracks", int(3))),
maxChi2ndf_(iConfig.getUntrackedParameter("maxChi2ndf", double(10.0)))
{
	produces<int>("flag");//0=invalid, 1=valid
	produces<reco::VertexCollection>("primVtx");//has to be vector. save one of length one
}


PrimVtxSelector::~PrimVtxSelector(){}

bool PrimVtxSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
	bool foundGoodVtx = false;
	iEvent_ = &iEvent;
	cnt++;

	std::auto_ptr<reco::VertexCollection> primaryVertex_ = std::auto_ptr<reco::VertexCollection>(new reco::VertexCollection);
	reco::VertexCollection & primaryVertex = * primaryVertex_;

	if(checkPrimVtx(primaryVertex)) foundGoodVtx = true;
	
	iEvent_->put(primaryVertex_,"primVtx");

	std::auto_ptr<int> flagPtr = std::auto_ptr<int>(new int);
	int &flag = *flagPtr;
	if(foundGoodVtx) flag = 1;
	else flag = 0;
	iEvent_->put(flagPtr,"flag");

	return foundGoodVtx;
}

// ------------ method called once each job just before starting event loop  ------------
void PrimVtxSelector::beginJob(){
	cnt = 0;
	cntFound = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void PrimVtxSelector::endJob() {
	float ratio = 0.0;
	if(cnt!=0) ratio=(float)cntFound/cnt;
    edm::LogInfo("PrimVtxSelector")<<"--> [PrimVtxSelector] asks for vertex with >= "<<minTracks_<<" tracks and chi2ndf <= "<<maxChi2ndf_<<". Efficiency: "<<cntFound<<"/"<<cnt<<" = "<<std::setprecision(4)<<ratio*100.0<<"%";
}

bool PrimVtxSelector::checkPrimVtx(reco::VertexCollection & primaryVertex){
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);
	
	std::vector<const reco::Vertex*> vtx;
	LogTrace("PrimVtxSelector")<<"evt "<<iEvent_->id().event()<<", lum "<<iEvent_->id().luminosityBlock()<<", run "<<iEvent_->id().run()<<" PrimVtxSelector::checkPrimVtx: no. of primVtx = "<<primVtxs->size();
	for(reco::VertexCollection::const_iterator v = primVtxs->begin(); v != primVtxs->end(); ++v){		
		if(v->tracksSize() >= minTracks_ && v->normalizedChi2() <= maxChi2ndf_) vtx.push_back(&*v);
		if(primVtxs->size() > 0) LogTrace("PrimVtxSelector")<<"evt "<<iEvent_->id().event()<<" PrimVtxSelector::checkPrimVtx: #trks "<<v->tracksSize()<<", chi2 "<<v->chi2()<<", ndf "<<v->ndof()<<", ("<<v->x()<<", "<<v->y()<<", "<<v->z()<<")+-("<<v->xError()<<", "<<v->yError()<<", "<<v->zError()<<")";
	}
	
	if(vtx.size()<1){
		LogTrace("PrimVtxSelector")<<"evt "<<iEvent_->id().event()<<" PrimVtxSelector::checkPrimVtx: No valid primary vertex found. Skip event.";
		return false;
	}
//	if(vtx.size()>1){
//		sort(vtx.begin(), vtx.end(), cmpNormalizedChi2<const reco::Vertex*>);
//		LogTrace("PrimVtxSelector")<<"evt "<<iEvent_->id().event()<<" PrimVtxSelector::checkPrimVtx: More than one ("<<vtx.size()<<") primary vertex found. Select best chi2ndf of "<<vtx.at(0)->normalizedChi2()<<".";
//	}
	cntFound++;
	primaryVertex.push_back(*(vtx.front()));
	return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(PrimVtxSelector);
