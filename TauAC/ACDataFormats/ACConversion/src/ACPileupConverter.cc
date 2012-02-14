#include "../interface/ACPileupConverter.h"

ACPileupConverter::ACPileupConverter(const PileupSummaryInfo & pileup) {
    num_PU_vertices_ = pileup.getPU_NumInteractions();
    zpositions_ = pileup.getPU_zpositions();
    sumpT_lowpT_ = pileup.getPU_sumpT_lowpT();
    sumpT_highpT_ = pileup.getPU_sumpT_highpT();
    ntrks_lowpT_ = pileup.getPU_ntrks_lowpT();
    ntrks_highpT_ = pileup.getPU_ntrks_highpT();
    bunchCrossing_ = pileup.getBunchCrossing();
    TrueNumInteractions_ = pileup.getTrueNumInteractions();
    instLumi_ = pileup.getPU_instLumi();
}
