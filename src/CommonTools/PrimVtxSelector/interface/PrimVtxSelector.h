// -*- C++ -*-
//
// Package:    PrimVtxSelector
// Class:      PrimVtxSelector
// 
/**\class PrimVtxSelector PrimVtxSelector.cc HiggsKinTau/PrimVtxSelector/src/PrimVtxSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Lars Perchalla
//         Created:  Thu Nov 12 14:10:26 CET 2009
// $Id: PrimVtxSelector.h,v 1.1 2010/01/21 09:56:08 perchall Exp $
//
//


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
	unsigned int verbosity_;
	
	unsigned int cnt, cntFound;
	
	template <typename T> static bool cmpNormalizedChi2(const T &a, const T &b){
		return a->normalizedChi2() > b->normalizedChi2();
	}
	
};
