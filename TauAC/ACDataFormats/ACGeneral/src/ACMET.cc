#include "../interface/ACMET.h"

ClassImp(ACMET)

ACMET::ACMET() {
    px_ = -1000.;
    py_ = -1000.;
    sumEt_ = -1000.;
    significanceMatrix_.ResizeTo(TMatrixD(2,2));
}
ACMET::ACMET(double px, double py, double sumEt, const TMatrixD & significanceMatrix) {
    px_ = px;
    py_ = py;
    sumEt_ = sumEt;
    significanceMatrix_.ResizeTo(TMatrixD(2,2));
    significanceMatrix_ = significanceMatrix;
}

/// the following functions are meant to not directly access the data members to allow for easy scaling of the values
const TVector3 ACMET::momentum() const { return TVector3(px(), py(), 0.); }
double ACMET::p() const { return Mag(); }
double ACMET::energy() const { return Mag();}
double ACMET::et() const { return Mag();}
double ACMET::pt() const { return Mag();}
double ACMET::phi() const { return momentum().Phi();}
double ACMET::theta() const { return momentum().Theta();}
double ACMET::eta() const { return momentum().Eta();}
const double ACMET::significance() const {
    bool filledMatrix = false;
    TMatrixD metmat = significanceMatrix();//need a non const copy here
    for (unsigned int i = 0; i!= 2; i++) {
        for (unsigned int j = 0; j!= 2; j++) {
            if (metmat(i,j)!=0.) {
                filledMatrix = true;
                continue;
            }
        }
    }
    if (!filledMatrix) return -1.;
    TVectorD metvec(2);
    metvec(0) = px();
    metvec(1) = py();
    double signif = -1.;
    if(std::fabs(metmat.Determinant())>0.000001){
        metmat.Invert();
        signif = metvec * (metmat * metvec);
    }
    return signif;
}