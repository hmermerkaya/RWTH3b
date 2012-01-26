#ifndef Ztotautau_ControlSample_h
#define Ztotautau_ControlSample_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Ztotautau_ControlSample : public Selection {

 public:
  Ztotautau_ControlSample(TString Name_, TString id_);
  virtual ~Ztotautau_ControlSample();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     hasTag,
	     TagPt,
	     TagIso,
	     JetPt,
	     deltaPhi,
	     MET,
	     MT,
	     ZMassV,
	     charge,
	     TauIsRef,
	     TauIsIso,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> PmuoverEtau;
  std::vector<TH1D> PmuoverEtau_hplus;
  std::vector<TH1D> PmuoverEtau_hminus;

};
#endif
