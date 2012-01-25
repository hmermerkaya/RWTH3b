#include "../interface/ACFitParticleConverter.h"

ACFitParticleConverter::ACFitParticleConverter(const SelectedKinematicParticle & selPart, int decayCharge, const ACGenParticleRef * genRef) {
    status_ = selPart.status();
    name_ = selPart.name();
    if (selPart.name()!="neutrino") pdgId_ = getPdgId(selPart.name(), selPart.charge());
    else pdgId_ = getPdgId(selPart.name(), decayCharge);//to destinguish between nu and antinu
    charge_ = selPart.charge();
    p4_ = selPart.p4();
    vertex_ = selPart.vertex();
    matrix_ = selPart.matrix();
    if (genRef) genRef_ = *genRef;

    //not filled yet:    decayIndex_ = ???;
    //members of SelectedKinematicParticle not used yet:        reco::PFCandidateRef PFCandRef() const;

    initialP4_.SetVectM(TVector3(selPart.input_parameters()[3],selPart.input_parameters()[4],selPart.input_parameters()[5]), selPart.input_parameters()[6]);
    initialVertex_.SetXYZ(selPart.input_parameters()[0],selPart.input_parameters()[1],selPart.input_parameters()[2]);
    initialMatrix_.ResizeTo(TMatrixDSym(7));
    initialMatrix_ = selPart.input_matrix();
}

int ACFitParticleConverter::getPdgId(const std::string name, const int charge) const {
    int pdgId = 0;
    if (name == "tau") pdgId = (int) charge * (-15);
    if (name == "pion") pdgId = (int) charge * (211);
    if (name == "neutrino") pdgId = (int) charge * (-16);//in case of neutrino the charge of the corresponding tau has to be provided

    return pdgId;
}