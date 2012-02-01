#include "TauSpinExample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

TauSpinExample::TauSpinExample(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TauSpinExample::~TauSpinExample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TauSpinExample::~TauSpinExample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TauSpinExample::~TauSpinExample()" << std::endl;
}

void  TauSpinExample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isZtautauto3pimu)          cut.at(isZtautauto3pimu)=1;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==isZtautauto3pimu){
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi\\pi\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi#pi#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  //mu
  mu_PmuoverEtau=HConfig.GetTH1D(Name+"_mu_PmuoverEtau","PmuoverEtau",20,0.0,1.0,"P_{#mu}/P_{#tau}","Events");
  mu_PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_mu_PmuoverEtau_hplus","PmuoverEtau_hplus",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{h^{+}}","Events");
  mu_PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_mu_PmuoverEtau_hminus","PmuoverEtau_hminus",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{h^{-}}","Events");
  mu_PmuoverEtau_Spin=HConfig.GetTH1D(Name+"_mu_PmuoverEtau_Spin","PmuoverEtau_Spin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{Spin}","Events");
  mu_PmuoverEtau_UnSpin=HConfig.GetTH1D(Name+"_mu_PmuoverEtau_UnSpin","PmuoverEtau_UnSpin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{UnSpin}","Events");
  mu_PmuoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_mu_PmuoverEtau_FlipSpin","PmuoverEtau_FlipSpin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{FlipSpin}","Events");
  mu_WT_Spin=HConfig.GetTH1D(Name+"_mu_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#mu}","Events");
  mu_WT_UnSpin=HConfig.GetTH1D(Name+"_mu_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#mu}","Events");
  mu_WT_FlipSpin=HConfig.GetTH1D(Name+"_mu_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#mu}","Events");

  //pinu
  pi_PmuoverEtau=HConfig.GetTH1D(Name+"_pi_PmuoverEtau","PmuoverEtau",20,0.0,1.0,"P_{#pi}/P_{#tau}","Events");
  pi_PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_PmuoverEtau_hplus","PmuoverEtau_hplus",20,0.0,1.0,"P_{#pi}/P_{#tau}|_{h^{+}}","Events");
  pi_PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_PmuoverEtau_hminus","PmuoverEtau_hminus",20,0.0,1.0,"P_{#pi}/P_{#tau}|_{h^{-}}","Events");
  pi_PmuoverEtau_Spin=HConfig.GetTH1D(Name+"_pi_PmuoverEtau_Spin","PmuoverEtau_Spin",20,0.0,1.0,"P_{#pi}/P_{#tau}|_{Spin}","Events");
  pi_PmuoverEtau_UnSpin=HConfig.GetTH1D(Name+"_pi_PmuoverEtau_UnSpin","PmuoverEtau_UnSpin",20,0.0,1.0,"P_{#pi}/P_{#tau}|_{UnSpin}","Events");
  pi_PmuoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_pi_PmuoverEtau_FlipSpin","PmuoverEtau_FlipSpin",20,0.0,1.0,"P_{#pi}/P_{#tau}|_{FlipSpin}","Events");
  pi_WT_Spin=HConfig.GetTH1D(Name+"_pi_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#pi}","Events");
  pi_WT_UnSpin=HConfig.GetTH1D(Name+"_pi_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#pi}","Events");
  pi_WT_FlipSpin=HConfig.GetTH1D(Name+"_pi_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#pi}","Events");

  //a1nu
  a1_PmuoverEtau=HConfig.GetTH1D(Name+"_a1_PmuoverEtau","PmuoverEtau",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}","Events");
  a1_PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_a1_PmuoverEtau_hplus","PmuoverEtau_hplus",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}|_{h^{+}}","Events");
  a1_PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_a1_PmuoverEtau_hminus","PmuoverEtau_hminus",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}|_{h^{-}}","Events");
  a1_PmuoverEtau_Spin=HConfig.GetTH1D(Name+"_a1_PmuoverEtau_Spin","PmuoverEtau_Spin",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}|_{Spin}","Events");
  a1_PmuoverEtau_UnSpin=HConfig.GetTH1D(Name+"_a1_PmuoverEtau_UnSpin","PmuoverEtau_UnSpin",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}|_{UnSpin}","Events");
  a1_PmuoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_a1_PmuoverEtau_FlipSpin","PmuoverEtau_FlipSpin",20,0.0,1.0,"P_{a_{1}(1260)}/P_{#tau}|_{FlipSpin}","Events");
  a1_WT_Spin=HConfig.GetTH1D(Name+"_a1_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{a_{1}(1260)}","Events");
  a1_WT_UnSpin=HConfig.GetTH1D(Name+"_a1_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{a_{1}(1260)}","Events");
  a1_WT_FlipSpin=HConfig.GetTH1D(Name+"_a1_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{a_{1}(1260)}","Events");

  //rhonu
  rho_PmuoverEtau=HConfig.GetTH1D(Name+"_rho_PmuoverEtau","PmuoverEtau",20,0.0,1.0,"P_{#rho}/P_{#tau}","Events");
  rho_PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_rho_PmuoverEtau_hplus","PmuoverEtau_hplus",20,0.0,1.0,"P_{#rho}/P_{#tau}|_{h^{+}}","Events");
  rho_PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_rho_PmuoverEtau_hminus","PmuoverEtau_hminus",20,0.0,1.0,"P_{#rho}/P_{#tau}|_{h^{-}}","Events");
  rho_PmuoverEtau_Spin=HConfig.GetTH1D(Name+"_rho_PmuoverEtau_Spin","PmuoverEtau_Spin",20,0.0,1.0,"P_{#rho}/P_{#tau}|_{Spin}","Events");
  rho_PmuoverEtau_UnSpin=HConfig.GetTH1D(Name+"_rho_PmuoverEtau_UnSpin","PmuoverEtau_UnSpin",20,0.0,1.0,"P_{#rho}/P_{#tau}|_{UnSpin}","Events");
  rho_PmuoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_rho_PmuoverEtau_FlipSpin","PmuoverEtau_FlipSpin",20,0.0,1.0,"P_{#rho}/P_{#tau}|_{FlipSpin}","Events");
  rho_WT_Spin=HConfig.GetTH1D(Name+"_rho_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#rho}","Events");
  rho_WT_UnSpin=HConfig.GetTH1D(Name+"_rho_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#rho}","Events");
  rho_WT_FlipSpin=HConfig.GetTH1D(Name+"_rho_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#rho}","Events");

  LongitudinalPolarization=HConfig.GetTH1D(Name+"_LongitudinalPolarization","LongitudinalPolarization",20,-2.0,2.0,"#rho_{l}","Events");
  LongitudinalPolarization_Spin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_Spin","LongitudinalPolarization_Spin",20,-2.0,2.0,"#rho_{l}|_{Spin}","Events");
  LongitudinalPolarization_UnSpin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_UnSpin","LongitudinalPolarization_UnSpin",20,-2.0,2.0,"#rho_{l}|_{UnSpin}","Events");
  LongitudinalPolarization_FlipSpin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_FlipSpin","LongitudinalPolarization_FlipSpin",20,-2.0,2.0,"#rho_{l}|_{FlipSpin}","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  TauSpinExample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&LongitudinalPolarization);
 Extradist1d.push_back(&LongitudinalPolarization_Spin);
 Extradist1d.push_back(&LongitudinalPolarization_UnSpin);
 Extradist1d.push_back(&LongitudinalPolarization_FlipSpin);

 Extradist1d.push_back(&mu_PmuoverEtau);
 Extradist1d.push_back(&mu_PmuoverEtau_hplus);
 Extradist1d.push_back(&mu_PmuoverEtau_hminus);
 Extradist1d.push_back(&mu_PmuoverEtau_Spin);
 Extradist1d.push_back(&mu_PmuoverEtau_UnSpin);
 Extradist1d.push_back(&mu_PmuoverEtau_FlipSpin);
 Extradist1d.push_back(&mu_WT_Spin);
 Extradist1d.push_back(&mu_WT_UnSpin);
 Extradist1d.push_back(&mu_WT_FlipSpin);

 Extradist1d.push_back(&pi_PmuoverEtau);
 Extradist1d.push_back(&pi_PmuoverEtau_hplus);
 Extradist1d.push_back(&pi_PmuoverEtau_hminus);
 Extradist1d.push_back(&pi_PmuoverEtau_Spin);
 Extradist1d.push_back(&pi_PmuoverEtau_UnSpin);
 Extradist1d.push_back(&pi_PmuoverEtau_FlipSpin);
 Extradist1d.push_back(&pi_WT_Spin);
 Extradist1d.push_back(&pi_WT_UnSpin);
 Extradist1d.push_back(&pi_WT_FlipSpin);

 Extradist1d.push_back(&a1_PmuoverEtau);
 Extradist1d.push_back(&a1_PmuoverEtau_hplus);
 Extradist1d.push_back(&a1_PmuoverEtau_hminus);
 Extradist1d.push_back(&a1_PmuoverEtau_Spin);
 Extradist1d.push_back(&a1_PmuoverEtau_UnSpin);
 Extradist1d.push_back(&a1_PmuoverEtau_FlipSpin);
 Extradist1d.push_back(&a1_WT_Spin);
 Extradist1d.push_back(&a1_WT_UnSpin);
 Extradist1d.push_back(&a1_WT_FlipSpin);

 Extradist1d.push_back(&rho_PmuoverEtau);
 Extradist1d.push_back(&rho_PmuoverEtau_hplus);
 Extradist1d.push_back(&rho_PmuoverEtau_hminus);
 Extradist1d.push_back(&rho_PmuoverEtau_Spin);
 Extradist1d.push_back(&rho_PmuoverEtau_UnSpin);
 Extradist1d.push_back(&rho_PmuoverEtau_FlipSpin);
 Extradist1d.push_back(&rho_WT_Spin);
 Extradist1d.push_back(&rho_WT_UnSpin);
 Extradist1d.push_back(&rho_WT_FlipSpin);


}

void  TauSpinExample::doEvent(){
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,tau_idx;
  value.at(isZtautauto3pimu)=0;
  pass.at(isZtautauto3pimu) = Ntp->NMCTaus()>=2;
  if(pass.at(isZtautauto3pimu))value.at(isZtautauto3pimu)=1;

  double wobs=1;
  double w=1;
  if(verbose)  std::cout << Ntp->GetMCID() << " " << Npassed.size() << " " << t << " " << Boson_idx << " " << " " << tau_idx 
	    << " " << Ntp->NMCTaus() << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);

    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);

    LongitudinalPolarization.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w);
    LongitudinalPolarization_Spin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*Spin_WT);
    LongitudinalPolarization_UnSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*UnSpin_WT);
    LongitudinalPolarization_FlipSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*FlipSpin_WT);
  
      ////////////////////////////////////////////////
      //
      // Spin Validation
      //
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_MUON,tau_idx)){
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
                ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	mu_PmuoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	mu_PmuoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	mu_PmuoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	mu_PmuoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	mu_PmuoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*Spin_WT);
	mu_PmuoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*Spin_WT);
	mu_WT_Spin.at(t).Fill(Spin_WT,w);
	mu_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	mu_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_PION,tau_idx)){
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
                ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	pi_PmuoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	pi_PmuoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	pi_PmuoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	pi_PmuoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	pi_PmuoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*Spin_WT);
	pi_PmuoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*Spin_WT);
	pi_WT_Spin.at(t).Fill(Spin_WT,w);
	pi_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	pi_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_A1_3PI,tau_idx)){
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
                ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	a1_PmuoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	a1_PmuoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	a1_PmuoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	a1_PmuoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	a1_PmuoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*Spin_WT);
	a1_PmuoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*Spin_WT);
	a1_WT_Spin.at(t).Fill(Spin_WT,w);
	a1_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	a1_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_RHO_PIPI0,tau_idx)){
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
	  if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
	    Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
	  else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_plus) ||
		  abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
		  abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
		  ){
	    X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
      }
      if(Tau_LV.E()>0){
	Tau_LV.Boost(-Boson_LV.BoostVector());
	X_LV.Boost(-Boson_LV.BoostVector());
	Boson_LV.Boost(-Boson_LV.BoostVector());
	// Now fill results                                                                                                                                                                                                                  
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	rho_PmuoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	rho_PmuoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	rho_PmuoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	rho_PmuoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	rho_PmuoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*Spin_WT);
	rho_PmuoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*Spin_WT);
	rho_WT_Spin.at(t).Fill(Spin_WT,w);
	rho_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	rho_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
  }
  if(verbose)std::cout << "done" << std::endl;
}




