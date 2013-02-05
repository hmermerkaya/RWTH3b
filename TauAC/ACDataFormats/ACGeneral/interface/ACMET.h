#ifndef ACMET_h
#define ACMET_h
/**
 Class to store missing Et objects.

 @author Lars Perchalla
 @date 2013
 */

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>

#include "ACCommon.h"

class ACMET : public ACCommon {
public:
    virtual std::string classname() const { return "ACMET"; }
    ACMET();
    ACMET(double px, double py, double sumEt, const TMatrixD & significanceMatrix);
    virtual ~ACMET() {};

    double px() const { return px_; }
    double py() const { return py_; }
    double sumEt() const { return sumEt_; }
    const TMatrixD & significanceMatrix() const { return significanceMatrix_; }

    /// get function compatible with ROOT (for convenience) 
    const double Mag() const { return sqrt(px_*px_ + py_*py_); }
    const TVector3 momentum() const;
    double p() const;
    double energy() const;
    double et() const;
    double pt() const;
    double phi() const;
    double theta() const;
    double eta() const;
    /// pseudo significance (NOT the real significance, but MET Significance = MET/sqrt(SumET))
    const double mEtSig() const { return ( sumEt_ ? (et()/sqrt(sumEt_)) : (0.0) ); }
    /// real MET significance
    const double significance() const;
    /// MET topology (topology = MET/SumET)
    const double topology() const { return ( sumEt_ ? (et()/sumEt_) : (0.0) ); }
    
    ClassDef(ACMET,1);

protected:
    /// momentum component x
    double px_;
    /// momentum component y
    double py_;
    /// total transverse energy
    double sumEt_;
    /// MET significance matrix 2x2 (may be unsymmetric, it's no covariance)
    TMatrixD significanceMatrix_;
};

typedef ACMET* pACMET;

#endif
