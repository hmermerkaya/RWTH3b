/**

 Simple primary vertex selector by minimum of included tracks and goodness of fit. produces one single vertex if possible.

 @author Lars Perchalla & Philip Sauerland
 @date 2010
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class PrimVtxSelector : public edm::EDFilter {
   public:
      explicit PrimVtxSelector(const edm::ParameterSet&);
      ~PrimVtxSelector();

private:
	virtual void beginJob() ;
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	bool checkPrimVtx(reco::VertexCollection & primaryVertex);
      
	edm::Event * iEvent_;
	edm::InputTag primVtx_;	
	unsigned int minTracks_;
	double maxChi2ndf_;
	
	unsigned int cnt, cntFound;
	
	template <typename T> static bool cmpNormalizedChi2(const T &a, const T &b){
		return a->normalizedChi2() > b->normalizedChi2();
	}
	
};
