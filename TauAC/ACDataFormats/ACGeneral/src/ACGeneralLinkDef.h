#include "../interface/ACCommon.h"
#include "../interface/ACTools.h"
#include "../interface/ACRef.h"
#include "../interface/ACCandidate.h"
#include "../interface/ACEventGlobals.h"
#include "../interface/ACEventInfo.h"
#include "../interface/ACGenParticle.h"
#include "../interface/ACJet.h"
#include "../interface/ACPFTau.h"
#include "../interface/ACParticle.h"
#include "../interface/ACFitParticle.h"
#include "../interface/ACVertex.h"
#include "../interface/ACDecay.h"
#include "../interface/ACGenDecay.h"
#include "../interface/ACFittedDecay.h"
#include "../interface/ACFittedThreeProngDecay.h"
#include "../interface/ACPileupInfo.h"
#include "../interface/ACTrigger.h"

#ifdef __CINT__

/// additional typedefs needed in the LinkDef file. all instances of different template types need to be linked.
typedef ACDecayBase<ACGenParticle> ACGenParticleDecayForLinkDefOnly;
typedef ACGenParticleDecayForLinkDefOnly* pACGenParticleDecayForLinkDefOnly;
typedef ACDecayBase<ACFitParticle> ACFitParticleDecayForLinkDefOnly;
typedef ACFitParticleDecayForLinkDefOnly* pACFitParticleDecayForLinkDefOnly;

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ACCommon+;
#pragma link C++ class ACTools+;

#pragma link C++ class ACCandidateRef+;
#pragma link C++ class ACGenParticleRef+;
#pragma link C++ class ACParticleRef+;
#pragma link C++ class ACFitParticleRef+;
#pragma link C++ class ACJetRef+;
#pragma link C++ class ACPFTauRef+;
#pragma link C++ class ACVertexRef+;
#pragma link C++ class ACGenDecayRef+;
#pragma link C++ class ACFittedDecayRef+;
#pragma link C++ class ACFittedThreeProngDecayRef+;

#pragma link C++ class ACCandidate+;
#pragma link C++ class pACCandidate+;
#pragma link C++ class std::vector<ACCandidate >;
#pragma link C++ class std::vector<pACCandidate >;

#pragma link C++ class ACPFTau+;
#pragma link C++ class pACPFTau+;
#pragma link C++ class std::vector<ACPFTau >;
#pragma link C++ class std::vector<pACPFTau >;

#pragma link C++ class ACEventGlobals+;
#pragma link C++ class pACEventGlobals+;

#pragma link C++ class ACEventInfo+;
#pragma link C++ class pACEventInfo+;

#pragma link C++ class ACGenParticle+;
#pragma link C++ class pACGenParticle+;
#pragma link C++ class std::vector<ACGenParticle >;
#pragma link C++ class std::vector<pACGenParticle >;
#pragma link C++ class ACParticle+;
#pragma link C++ class pACParticle+;
#pragma link C++ class std::vector<ACParticle >;
#pragma link C++ class std::vector<pACParticle >;
#pragma link C++ class ACFitParticle+;
#pragma link C++ class pACFitParticle+;
#pragma link C++ class std::vector<ACFitParticle >;
#pragma link C++ class std::vector<pACFitParticle >;

#pragma link C++ class std::pair<std::string, bool>+;

#pragma link C++ class ACGenParticleDecayForLinkDefOnly+;
#pragma link C++ class pACGenParticleDecayForLinkDefOnly+;
#pragma link C++ class ACGenDecay+;
#pragma link C++ class pACGenDecay+;
#pragma link C++ class std::vector<ACGenDecay >;
#pragma link C++ class std::vector<pACGenDecay >;
#pragma link C++ class ACDecay+;
#pragma link C++ class pACDecay+;
#pragma link C++ class std::vector<ACDecay >;
#pragma link C++ class std::vector<pACDecay >;
#pragma link C++ class ACFitParticleDecayForLinkDefOnly+;
#pragma link C++ class pACFitParticleDecayForLinkDefOnly+;
#pragma link C++ class ACFittedDecay+;
#pragma link C++ class pACFittedDecay+;
#pragma link C++ class std::vector<ACFittedDecay >;
#pragma link C++ class std::vector<pACFittedDecay >;
#pragma link C++ class ACFittedThreeProngDecay+;
#pragma link C++ class pACFittedThreeProngDecay+;
#pragma link C++ class std::vector<ACFittedThreeProngDecay >;
#pragma link C++ class std::vector<pACFittedThreeProngDecay >;

#pragma link C++ class ACJet+;
#pragma link C++ class pACJet+;
#pragma link C++ class std::vector<ACJet >;
#pragma link C++ class std::vector<pACJet >;

#pragma link C++ class ACVertex+;
#pragma link C++ class pACVertex+;
#pragma link C++ class std::vector<ACVertex >;
#pragma link C++ class std::vector<pACVertex >;

#pragma link C++ class ACPileupInfo+;

#pragma link C++ class ACTrigger+;
#pragma link C++ class pACTrigger+;

#endif
