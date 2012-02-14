/**
 CMSSW producer to filter different candidate collections.
 
 @author Lars Perchalla & Philip Sauerland
 @date 2012
 */


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"



class MultiCandidateSelector : public edm::EDProducer {
public:
    explicit MultiCandidateSelector(const edm::ParameterSet&);
    ~MultiCandidateSelector();
    
    /// method fills 'descriptions' with the allowed parameters for the module
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    /// collection of filters
    std::vector<edm::ParameterSet> filters_;
    
    /// register products
    void registerProducts();
    
    /// method called once each job just before starting event loop
    virtual void beginJob() ;
    /// method called once for each event
    virtual void produce(edm::Event&, const edm::EventSetup&);
    /// method called once each job just after ending the event loop
    virtual void endJob() ;
    /// method called when starting to processes a run
    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    /// method called when ending the processing of a run
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    /// method called when starting to processes a luminosity block
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    /// method called when ending the processing of a luminosity block
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
};
