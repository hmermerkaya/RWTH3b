#include "../interface/ACExampleAnalyzer.h"

ACExampleAnalyzer::ACExampleAnalyzer() {
}

int main(int argc, char* argv[]) {
    ACExampleAnalyzer anlzr;
    std::vector<ACTreeReader* > loops;
    std::vector<std::string> filenames;
    filenames.push_back("tauSMBrFromZ_-1_7TeV.root");
    loops.push_back(new ACTreeReader(filenames, "FinalTreeFiller/TauACEvent"));
    loops.back()->loop(anlzr, -1);

    return 0;
}

void ACExampleAnalyzer::analyze(const ACEvent & event) {
    printf("Working on event %d, run %d, and lumi %d. Type is '%s'\n", event.eventInfo()->eventID(), event.eventInfo()->runID(), event.eventInfo()->lumiID(), event.eventInfo()->type().c_str());

    printf("PFMET topology: %f\n", event.eventGlobals()->pfMETTopology());

    printf("trigger: %s\n", event.trigger()->tableName().c_str());
    printf("offlinePV: size %lu\n", event.offlinePV().size());
    printf("reducedPV: size %lu\n", event.reducedPV().size());

    printf("number of generator particles: %lu\n", event.generator().size());
    for (std::vector<ACGenParticle *>::const_iterator ip = event.generator().begin(); ip != event.generator().end(); ++ip) {
        printf("\t gen particle name=%s, pdgID=%d, pt %f\n", (*ip)->name().c_str(), (*ip)->pdgId(), (*ip)->pt());
        if ((*ip)->genDecayRef().isValid()) printf("\t gen tau decay ref name=%s\n", (*ip)->genDecayRef()->name().c_str());
        else printf("\t no gen tau decay ref\n");
    }
    printf("number of generator tau decays: %lu\n", event.genTauDecays().size());
    for (std::vector<ACGenDecay *>::const_iterator decay = event.genTauDecays().begin(); decay != event.genTauDecays().end(); ++decay) {
        printf("gen tau decay name=%s\n", (*decay)->name().c_str());
        if ((*decay)->particles()->size()>0) {
            printf("gen decay stores %lu particles.\n", (*decay)->particles()->size());
            for (std::vector<ACGenParticleRef>::const_iterator ip = (*decay)->particles()->begin(); ip != (*decay)->particles()->end(); ++ip) {
                printf("\t pdgID %i, pt %f\n", (*ip)->pdgId(), (*ip)->pt());
            }
        }
    }
    printf("number of muons: %lu\n", event.muons().size());
    for (std::vector<ACParticle *>::const_iterator ip = event.muons().begin(); ip != event.muons().end(); ++ip) {
        printf(" pt %f\n", (*ip)->pt());
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
        printf("\t\t PFTauRef pt: %f\n", (*decay)->PFTauRef()->p4().Pt());
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
        printf("\t tau pt: %f, charged hadr signal cands: %i\n", (*itau)->p4().Pt(), (*itau)->signalPFChargedHadrCands());
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
