// -*- C++ -*-
//
// Package:    TutorialAnlzr
// Class:      TutorialAnlzr
// 
/**
 *
 * Description:
 *
 * Implementation:
 * <Notes on implementation>
 */
//
// Original Author:  Lars Perchalla
//         Created:  Thu Jun 10 10:45:54 CEST 2010
// $Id: TutorialAnlzr.cc,v 1.3 2010/08/13 10:47:07 perchall Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>


class TutorialAnlzr : public edm::EDAnalyzer {
public:
	explicit TutorialAnlzr(const edm::ParameterSet&);
	~TutorialAnlzr();
	
	
private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	edm::InputTag generatorTag_;
};
