// -*- C++ -*-
//
// Package:    KinematicTau
// Class:      KinematicTau
// 
/**\class KinematicTau KinematicTau.cc RecoTauTag/KinematicTau/src/KinematicTau.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Philip Sauerland
//         Created:  Tue Jan 12 15:13:30 CET 2010
// $Id$
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
//
// class decleration
//

class KinematicTau : public edm::EDAnalyzer {
   public:
      explicit KinematicTau(const edm::ParameterSet&);
      ~KinematicTau();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KinematicTau::KinematicTau(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


KinematicTau::~KinematicTau()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
KinematicTau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
KinematicTau::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KinematicTau::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicTau);
