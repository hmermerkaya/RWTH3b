#include "../interface/PrimVtxSelector.h"

PrimVtxSelector::PrimVtxSelector(const edm::ParameterSet& iConfig):
primVtx_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),//primVtx from generalTracks
minTracks_(iConfig.getUntrackedParameter("minTracks", int(3))),
maxChi2ndf_(iConfig.getUntrackedParameter("maxChi2ndf", double(10.0))),
verbosity_(iConfig.getUntrackedParameter("verbosity", int(2)))
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
	printf("--> [PrimVtxSelector] asks for vertex with >= %i tracks and chi2ndf <= %f. Efficiency: %d/%d = %.2f%%\n", minTracks_, maxChi2ndf_, cntFound, cnt, ratio*100.0);
}

bool PrimVtxSelector::checkPrimVtx(reco::VertexCollection & primaryVertex){
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent_->getByLabel( primVtx_, primVtxs);
	
	std::vector<const reco::Vertex*> vtx;
	for(reco::VertexCollection::const_iterator v = primVtxs->begin(); v != primVtxs->end(); ++v){
		//printf("primVtx: chi2=%f, ndof=%f", v->chi2(), v->ndof());
		if(v->tracksSize() >= minTracks_ && v->normalizedChi2() <= maxChi2ndf_) vtx.push_back(&*v);
		if(verbosity_>=2) if(primVtxs->size() > 1) printf("evt %d PrimVtxSelector::checkPrimVtx: #trks %3d, chi2 %9.6f, ndf %5.3f, (%8.6f, %8.6f, %8.6f)+-(%8.6f, %8.6f, %8.6f)\n", iEvent_->id().event(), v->tracksSize(), v->chi2(), v->ndof(), v->x(), v->y(), v->z(), v->xError(), v->yError(), v->zError() );
	}
	
	if(vtx.size()<1){
		LogTrace("PrimVtxSelector")<<"PrimVtxSelector::checkPrimVtx: No valid primary vertex found. Skip event "<<iEvent_->id().event()<<".";
		return false;
	}
	if(vtx.size()>1){
		sort(vtx.begin(), vtx.end(), cmpNormalizedChi2<const reco::Vertex*>);
		if(verbosity_>=1) printf("evt %d PrimVtxSelector::checkPrimVtx: More than one (%i) primary vertex found. Select best chi2ndf of %f.\n", iEvent_->id().event(), vtx.size(), vtx.at(0)->normalizedChi2());
	}
	cntFound++;
	primaryVertex.push_back(*(vtx.front()));
	return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(PrimVtxSelector);
