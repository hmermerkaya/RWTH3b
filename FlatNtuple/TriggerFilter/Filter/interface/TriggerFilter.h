// -*- C++ -*-
//
// Package:    Filter
// Class:      Filter
// 
/**\class Filter Filter.cc TriggerFilter/Filter/src/Filter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Cherepanov
//         Created:  Wed Nov 23 16:24:12 CET 2011
// $Id: Filter.cc,v 1.2 2012/02/23 15:24:29 cherepan Exp $
//
//
#ifndef TriggerFilter_h
#define TriggerFilter_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TriggerFilter/Filter/interface/TriggerHelper.h"
//
// class declaration
//

class TriggerFilter : public edm::EDFilter {
  public:
    explicit TriggerFilter(const edm::ParameterSet&);
    ~TriggerFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
      
    virtual bool beginRun(edm::Run&, edm::EventSetup const&);
    virtual bool endRun(edm::Run&, edm::EventSetup const&);
    virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    //Filter Functions
    virtual bool FilteronTauplusXTrigger(edm::Event&);
    virtual bool FilteronMuonTrigger(edm::Event&);
    virtual bool FilteronElectronTrigger(edm::Event&);
    virtual bool FilteronTauplusMETTrigger(edm::Event&);

    // ----------member data ---------------------------
    edm::Event * iEvent_;
    std::string processName_; // process name of (HLT) process for which to get HLT configuration
    HLTConfigProvider hltConfig_;/// The instance of the HLTConfigProvider as a data member
    int cnt_;
    int cntFound_;
    TriggerHelper MyTriggerHelper;
    // FilterFlags
    bool doTauplusXTrigger_;
    bool doMuonTrigger_;
    bool doElectronTrigger_;
    bool doTauplusMETTrigger_;
};

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerFilter);
#endif
