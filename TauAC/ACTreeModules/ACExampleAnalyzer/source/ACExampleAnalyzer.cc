#include "../interface/ACExampleAnalyzer.h"

ACExampleAnalyzer::ACExampleAnalyzer() {
}

int main(int argc, char* argv[]) {
    
    /// one will want to parse a config file here to define steering parameters
    
    ACExampleAnalyzer anlzr;
    std::vector<ACTreeReader* > loops;
    std::vector<std::string> filenames;
    //filenames.push_back("/user/perchalla/output/analysis/CMSSW_4_4_2/development/SingleMuNov2011A_100_7TeV.root");
    filenames.push_back("/user/perchalla/output/analysis/CMSSW_4_4_2/development/tauSMBrFromZ_100_7TeV.root");    
    ACTreeReader::SetVerbosity(true);
    loops.push_back(new ACTreeReader(filenames, "FinalTreeFiller/TauACEvent"));
    loops.back()->loop(anlzr, -1);

    return 0;
}

void ACExampleAnalyzer::analyze(const ACEvent & event) {
    printf("Working on event %d, run %d, and lumi %d. Type is '%s'\n", event.eventInfo()->eventID(), event.eventInfo()->runID(), event.eventInfo()->lumiID(), event.eventInfo()->type().c_str());
    printf("Filter results:");
    for (std::vector<std::pair<std::string, int> >::const_iterator filter = event.eventInfo()->filterResults()->begin(); filter != event.eventInfo()->filterResults()->end(); ++filter) {
        printf("\t%s=%i", filter->first.c_str(), filter->second);
    }
    printf("\n");
    printf("PFMET topology: %f\n", event.eventGlobals()->pfMET().topology());
    printf("PFMET MET, METType1Corrected: %f, %f\n", event.eventGlobals()->pfMET().Mag(), event.eventGlobals()->pfType1CorrectedMET().Mag());

    const std::map<std::string, std::vector<std::string> > & moduleNamesPerPath = event.trigger()->moduleNamesPerPath();
    const std::map<std::string, std::vector<ACTriggerObjectRef> > & objectsPerModule = event.trigger()->objectsPerModule();
    printf("trigger: %s, with %lu paths and %lu trigger objects.\n", event.trigger()->tableName().c_str(), moduleNamesPerPath.size(), event.triggerObjects().size());
    for (std::map<std::string, std::vector<std::string> >::const_iterator path = moduleNamesPerPath.begin(); path != moduleNamesPerPath.end(); ++path) {
        unsigned int objectSize = 0;
        for (std::vector<std::string>::const_iterator module = path->second.begin(); module !=path->second.end(); ++module) {
            std::map<std::string, std::vector<ACTriggerObjectRef> >::const_iterator objects = objectsPerModule.find(*module);
            if (objects == objectsPerModule.end()) continue;
            objectSize += objects->second.size();
            if (objects->second.size() > 0) {
                printf("\t\tlast object of module %s: %i pt %f\n", module->c_str(), objects->second.back()->pdgID(), objects->second.back()->pt());
            }
        }
        printf("\tpath %s contains %lu modules and %i trigger objects\n", path->first.c_str(), path->second.size(), objectSize);
    }
    
    printf("beamspot:\n");
    std::stringstream output;
    event.beamSpot()->print(output);
    std::cout<<output.str();
    printf("offlinePV: size %lu\n", event.offlinePV().size());
    for (std::vector<ACVertex*>::const_iterator vtx = event.offlinePV().begin(); vtx != event.offlinePV().end(); ++vtx) {
        std::cout<<"\tstored "<<(*vtx)->tracks().size()<<"/"<<(*vtx)->trackSize();
        if ((*vtx)->tracks().size() > 0) std::cout<<" (first track with pt "<<(*vtx)->tracks().front()->pt()<<")";
        std::cout<<std::endl;
    }
    printf("reducedPV: size %lu\n", event.reducedPV().size());

    printf("number of generator particles: %lu\n", event.generator().size());
    for (std::vector<ACGenParticle *>::const_iterator ip = event.generator().begin(); ip != event.generator().end(); ++ip) {
        printf("\t gen particle name=%s, pdgID=%d, status=%i, pt %f", (*ip)->name().c_str(), (*ip)->pdgId(), (*ip)->status(), (*ip)->pt());
        if ((*ip)->mother().isValid()) printf("\t mother: gen particle name=%s, pdgID %i, pt %f", (*ip)->mother()->name().c_str(), (*ip)->mother()->pdgId(), (*ip)->mother()->pt());
        printf("\t daughters: %lu\n", (*ip)->daughters().size());
        if ((*ip)->genDecayRef().isValid()) printf("\t gen tau decay ref name=%s\n", (*ip)->genDecayRef()->name().c_str());
    }
    printf("number of generator tau decays: %lu\n", event.genTauDecays().size());
    for (std::vector<ACGenDecay *>::const_iterator decay = event.genTauDecays().begin(); decay != event.genTauDecays().end(); ++decay) {
        printf("gen tau decay name=%s\n", (*decay)->name().c_str());
        if ((*decay)->particles()->size()>0) {
            printf("gen decay stores %lu particles. There are %i stable and %i unstable daughters.\n", (*decay)->particles()->size(), (*decay)->numberOfDaughters(), (*decay)->numberOfUnstableDaughters());
            for (std::vector<ACGenParticleRef>::const_iterator ip = (*decay)->particles()->begin(); ip != (*decay)->particles()->end(); ++ip) {
                printf("\t pdgID %i, pt %f\n", (*ip)->pdgId(), (*ip)->pt());
            }
        }
    }
    printf("number of tracks: %lu\n", event.tracks().size());
    for (std::vector<ACTrack *>::const_iterator ip = event.tracks().begin(); ip != event.tracks().end(); ++ip) {
        printf(" pt %f\n", (*ip)->pt());
    }
    printf("number of muons: %lu\n", event.muons().size());
    for (std::vector<ACMuon *>::const_iterator ip = event.muons().begin(); ip != event.muons().end(); ++ip) {
        std::cout<<" pt "<<(*ip)->pt()<<", isGlobalMuonPromptTight "<<(*ip)->isGlobalMuonPromptTight()<<", isRecommendedMuon "<<(*ip)->isRecommendedMuon()<<", isolation "<<(*ip)->isolation()<<" (isIsolatedLoose "<<(*ip)->isIsolatedLoose()<<", isIsolatedTight "<<(*ip)->isIsolatedTight()<<")";
        if ((*ip)->trackRef().isValid()) std::cout<<"\ttrackRef "<<(*ip)->trackRef()->pt()<<std::endl;
        else std::cout<<"\tinvalid trackRef"<<std::endl;
    }
    printf("number of electrons: %lu\n", event.electrons().size());
    for (std::vector<ACParticle *>::const_iterator ip = event.electrons().begin(); ip != event.electrons().end(); ++ip) {
        printf(" pt %f\n", (*ip)->pt());
    }
    printf("number of fitted tau decays: %lu\n", event.tauDecays().size());
    for (std::vector<ACFittedThreeProngDecay *>::const_iterator decay = event.tauDecays().begin(); decay != event.tauDecays().end(); ++decay) {
        printf("decay stores %lu particles.\n", (*decay)->particles()->size());
        printf("\t\t PFTauRef pt: %f\n", (*decay)->PFTauRef()->p4().Pt());
        for (std::map<std::string,bool>::const_iterator idiscr = (*decay)->PFTauRef()->discriminators().begin(); idiscr != (*decay)->PFTauRef()->discriminators().end(); ++idiscr) {
            printf("\t\t\t discr: %s pass: %i\n", idiscr->first.c_str(), idiscr->second);
        }
        printf("\t SV (%f, %f, %f) error matrix:\n", (*decay)->tau()->vertex().X(), (*decay)->tau()->vertex().Y(), (*decay)->tau()->vertex().Z());
        (*decay)->tau()->matrix().Print();
        for (std::vector<ACFitParticleRef>::const_iterator ip = (*decay)->particles()->begin(); ip != (*decay)->particles()->end(); ++ip) {
            printf("\t pdgID %i, pt %f\n", (*ip)->pdgId(), (*ip)->pt());
            if ((*ip)->genRef().isValid()) {
                printf("\t\t ref to pdgID %i, pt %f\n", (*ip)->genRef()->pdgId(), (*ip)->genRef()->pt());
            }
        }
    }
    printf("number of pf jets: %lu\n", event.pfJets().size());
    for (std::vector<ACJet *>::const_iterator ijet = event.pfJets().begin(); ijet != event.pfJets().end(); ++ijet) {
        printf("\t jet energy: %f, charge: %f, and multiplicity: %i\n", (*ijet)->p4().E(), (*ijet)->charge(), (*ijet)->multiplicity());
    }
    printf("number of pf taus: %lu\n", event.pfTaus().size());
    for (std::vector<ACPFTau *>::const_iterator itau = event.pfTaus().begin(); itau != event.pfTaus().end(); ++itau) {
        if ((*itau)->jetRef().isValid()) {
            printf("\t tau pt: %f, charged hadr signal cands: %i, pi0 signal cands: %i, jetRef pt: %f\n", (*itau)->p4().Pt(), (*itau)->signalPFChargedHadrCands(), (*itau)->signalPiZeroCands(), (*itau)->jetRef()->pt());
        } else {
            printf("\t tau pt: %f, charged hadr signal cands: %i, pi0 signal cands: %i, no jetRef!\n", (*itau)->p4().Pt(), (*itau)->signalPFChargedHadrCands(), (*itau)->signalPiZeroCands());
        }
    }
    printf("number of pileup information: %lu\n", event.pileup().size());
    for (std::vector<ACPileupInfo *>::const_iterator pileup = event.pileup().begin(); pileup != event.pileup().end(); ++pileup) {
        printf("\t pileup: bunch crossing %i, number of pileup interactions %i\n", (*pileup)->getBunchCrossing(), (*pileup)->getPU_NumInteractions());
    }
    
    printf("\n");
}
void ACExampleAnalyzer::beginRun() {
    printf("This is the beginRun!\n");
}
void ACExampleAnalyzer::endRun() {
    printf("This is the endRun!\n");
}
void ACExampleAnalyzer::beginLuminosityBlock() {
    printf("This is the beginLuminosityBlock!\n");
}
void ACExampleAnalyzer::endLuminosityBlock() {
    printf("This is the endLuminosityBlock!\n");
}
void ACExampleAnalyzer::beginJob() {
    printf("This is the beginJob!\n");
}
void ACExampleAnalyzer::endJob() {
    printf("This is the endJob!\n");
}
