#include "Selection.h"

#include "Tables.h"
#include "Plots.h"
#include "SkimConfig.h"
#include <cstdlib>
#include <map>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

Selection::Selection(TString Name_, TString id_):
  Selection_Base(Name_,id_)
  ,NGoodFiles(0)
  ,NBadFiles(0)
  ,isStored(false)
  ,data(0)
  ,HConfig()
  ,Lumi(1)
{
  if(Name_)
  ListofBadFiles.clear();
}

Selection::~Selection(){
  std::cout << Get_Name() << " NGoodFile= " <<  NGoodFiles << " NBadFiles=" << NBadFiles << std::endl;
  if(ListofBadFiles.size()>0)std::cout <<  " List of Bad Files:" << std::endl;
  for(int i=0; i< ListofBadFiles.size(); i++){
    std::cout <<  ListofBadFiles.at(i) << std::endl;
  }
}

void Selection::ConfigureHistograms(){
  if(!isStored){
    if(verbose) std::cout << "Selection::ConfigureHistograms() starting" << std::endl;
    isStored=true;
    Store_ExtraDist();
    for(unsigned int i=0; i<Nminus1.size(); i++){
      for(unsigned int j=0; j<Nminus1.at(i).size();j++){
	if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point A j= " << j << std::endl;
	Nminus1.at(i).at(j).Sumw2();
	Nminus0.at(i).at(j).Sumw2();
	if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point B j= " << j << std::endl;
	if(distindx.at(i)){
	  Nminus1dist.at(i).at(j).Sumw2();
	  Accumdist.at(i).at(j).Sumw2();
	}
	if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point C j= " << j << std::endl;
	if(i==0){
	  if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point 1 j= " << j << std::endl;
	  TString name=Npassed.at(j).GetName();
	  name+="_noweight";
	  Npassed_noweight.push_back((*((TH1D*)Npassed.at(j).Clone(name))));
	  Npassed_noweight.at(j).Sumw2();
	  Npassed.at(j).Sumw2();
	  if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point 2 j= " << j << std::endl;
	  for(unsigned int k=0; k<Extradist1d.size();k++){
	    if(verbose) std::cout << "k 1d" << k <<  std::endl;
	    if(verbose) std::cout << "k " << k <<  " " << j << " " << Extradist1d.at(k)->at(j).GetName() << std::endl;
	    Extradist1d.at(k)->at(j).Sumw2();
	  }
	  if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point 3 j= " << j
				<< " size= " << Extradist2d.size() << std::endl;
	  for(unsigned int k=0; k<Extradist2d.size();k++){
	    if(verbose) std::cout << "k 2d" << k <<  std::endl;
	    if(verbose) std::cout << "k " << k <<  " " << j << " " << Extradist2d.at(k)->size() << " " 
				  <<  Extradist2d.at(k)->at(j).GetName() << std::endl;
	    Extradist2d.at(k)->at(j).Sumw2();
	    if(verbose) std::cout << "OK " << std::endl;
	  }
	  if(verbose) std::cout << "Selection::ConfigureHistograms() i= " << i << "  Point 4 j= " << j << std::endl;
	}
      }
    }
    if(verbose) std::cout << "Selection::ConfigureHistograms() done" << std::endl;
  }
}

void Selection::LoadResults(std::vector<TString> files){ 
  if(!isStored){
    ConfigureHistograms();
  }
  for(int f=0;f<files.size();f++){
    TString file=files.at(f);
    if(!file.Contains(".root")){
      vector<TString> filelist;
      string dir=file.Data();
      DIR *dp;
      struct dirent *dirp;
      if((dp  = opendir(dir.c_str())) == NULL) {
	cout << "Error(" << errno << ") opening " << dir << endl;
      }
      else{
	while ((dirp = readdir(dp)) != NULL) {
	  filelist.push_back(string(dirp->d_name));
	}
	closedir(dp);
      }
      TString ID=Get_Name()+".root";
      for(int i=0;i<filelist.size();i++){
	if(filelist.at(i).Contains(ID)){
	  file+=filelist.at(i);
	  break;
	}
      }
    }
    TFile *f=TFile::Open(file,"READ");
    std::cout << "Selection::LoadResults " << file << std::endl;
    TString hname;
    if(f->IsOpen()){
      for(unsigned int i=0; i<Nminus1.size(); i++){
	for(unsigned int j=0; j<Nminus1.at(i).size();j++){
	  hname=(Nminus1.at(i).at(j)).GetName();
	  Nminus1.at(i).at(j).Add((TH1*)f->Get(hname),1.000);
	  hname=(Nminus0.at(i).at(j)).GetName();
	  Nminus0.at(i).at(j).Add((TH1*)f->Get(hname),1.000);
	  if(distindx.at(i)){
	    hname=(Nminus1dist.at(i).at(j)).GetName();
	    Nminus1dist.at(i).at(j).Add((TH1*)f->Get(hname),1.000);
	    hname=(Accumdist.at(i).at(j)).GetName();
	    Accumdist.at(i).at(j).Add((TH1*)f->Get(hname),1.000);
	  }
	  if(i==0){
	    hname=((Npassed.at(j)).GetName());
	    TH1* temp=(TH1*)f->Get(hname);
	    Npassed.at(j).Add(temp,1.000);
	    hname=((Npassed_noweight.at(j)).GetName());
	    TH1* tempnw=(TH1*)f->Get(hname);
	    Npassed_noweight.at(j).Add(tempnw,1.000);
	    for(unsigned int k=0; k<Extradist1d.size();k++){
	      hname=(Extradist1d.at(k)->at(j)).GetName();
	      Extradist1d.at(k)->at(j).Add((TH1*)f->Get(hname),1.000);
	    }
	    for(unsigned int k=0; k<Extradist2d.size();k++){
	      TString n=Extradist2d.at(k)->at(j).GetName();
	      if(!n.Contains("egammaMap")){
		hname=(Extradist2d.at(k)->at(j)).GetName();
		Extradist2d.at(k)->at(j).Add((TH1*)f->Get(hname),1.000);
	      }
	    }
	  }
	}
      }
      NGoodFiles++;
    }
    else{
      NBadFiles++;
      std::cout << "WARNING: " << file << " NOT OPENED" << std::endl;
      ListofBadFiles.push_back(file);
    }
    f->Close();
  }
}

bool Selection::AnalysisCuts(int t,double w,double wobjs){
  int ncuts=Nminus1.size();
  if(Npassed.size()!=Npassed_noweight.size()){
    std::cout << "ERROR Histograms not Configured. Please fix your code!!!! Running Selection::ConfigureHistograms()" << std::endl;
    Selection::ConfigureHistograms();
  }
   if(0<=t && t<types.size()){
    int nfail=0;
    int fail=-1;
    Npassed.at(t).Fill(-0.5,w);
    Npassed_noweight.at(t).Fill(-0.5,1);
    for(int i=0; i<ncuts;i++){
      if(!pass.at(i)){
	fail=i;
	nfail++;
      }
      if(nfail==0){
	Npassed.at(t).Fill((float)i+0.5,w*wobjs);
	Npassed_noweight.at(t).Fill((float)i+0.5,1);
	if(i+1<ncuts){
	  if(distindx.at(i+1)){
	    for(int k=0;k<dist.at(i+1).size();k++){
	      Accumdist.at(i+1).at(t).Fill(dist.at(i+1).at(k),w*wobjs);
	    }
	  }
	}
	if(i==0){
	  if(distindx.at(i)){
	    for(int k=0;k<dist.at(i).size();k++){
	      Accumdist.at(i).at(t).Fill(dist.at(i).at(k),w*wobjs);
	    }
	  }
	}
      }
    }
    if(nfail<=1){
      for(int i=0; i<ncuts;i++){
	if(fail==i || nfail==0){
	  Nminus1.at(i).at(t).Fill(value.at(i),w*wobjs);
	  if(distindx.at(i)){
	    for(int k=0;k<dist.at(i).size();k++){
	      Nminus1dist.at(i).at(t).Fill(dist.at(i).at(k),w*wobjs);
	    }
	  }
	}
      }
      
      if(nfail==0){
	for(int i=0; i<ncuts;i++){
	  Nminus0.at(i).at(t).Fill(value.at(i),w*wobjs);
	}
	return true;
      }
    }
  }
   std::cout << "Failed " << std::endl;
  return false;
}

void  Selection::Finish(){
  if(Npassed.size()!=Npassed_noweight.size()){
    std::cout << "ERROR Histograms not Configured. Please fix your code!!!! Running Selection::ConfigureHistograms()" << std::endl;
    Selection::ConfigureHistograms();
  }
  if(!isStored){
    ConfigureHistograms();
  }
  std::cout << "Writing out "+Name+".root ..." << std::endl;
  TString fName;
  if(runtype==GRID)         fName="GRID_";
  if(runtype==Local)        fName="LOCAL_";
  if(mode==RECONSTRUCT)     fName+="COMBINED_";
  if(mode==ANALYSIS)        fName+="ANALYSIS_";
  fName+=Name;

  TFile f(fName+".root","RECREATE");
  for(unsigned int i=0; i<Nminus1.size(); i++){
    for(unsigned int j=0; j<Nminus1.at(i).size();j++){
      Nminus1.at(i).at(j).Write((Nminus1.at(i).at(j)).GetName());
      Nminus0.at(i).at(j).Write((Nminus0.at(i).at(j)).GetName());
      if(distindx.at(i)){
	Nminus1dist.at(i).at(j).Write((Nminus1dist.at(i).at(j)).GetName());
	Accumdist.at(i).at(j).Write((Accumdist.at(i).at(j)).GetName());
      }
      if(i==0){
	Npassed.at(j).Write((Npassed.at(j)).GetName());
	Npassed_noweight.at(j).Write((Npassed_noweight.at(j)).GetName());
	for(unsigned int i=0; i<Extradist1d.size();i++){
	  Extradist1d.at(i)->at(j).Write((Extradist1d.at(i)->at(j)).GetName());
	}
	for(unsigned int i=0; i<Extradist2d.size();i++){
	  Extradist2d.at(i)->at(j).Write((Extradist2d.at(i)->at(j)).GetName());
	}
      }
    }
  }
  f.Close();
  std::vector<float> nevents;
  std::vector<float> nevents_noweight;

  SkimConfig SC;
  //SC.Load();
  //TString SkimFile=Name;
  //SC.LoadSkimEff("Tools/Zee_MET_CS_0SkimEff.dat");
  //SC.ApplySkimEfficiency(types,Npassed,Npassed_noweight);
  std::cout << "F" << std::endl;
  for(int i=0; i<Npassed.size();i++){
    nevents.push_back(Npassed.at(i).GetBinContent(1));
    nevents_noweight.push_back(Npassed_noweight.at(i).GetBinContent(1));
    std::cout << "Weights: " << Npassed.at(i).GetBinContent(1) << std::endl;
  }
  std::cout << "G" << std::endl;
  if(runtype!=GRID){
    std::cout << "Printing Plots " << std::endl;
    
    Plots P;
    P.Plot1D(Nminus1,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
    P.Plot1D(Nminus0,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
    P.Plot1D(Nminus1dist,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
    P.Plot1D(Accumdist,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
    
    for(unsigned int i=0; i<Extradist1d.size();i++){
      P.Plot1D((*Extradist1d.at(i)),Lumi,CrossSectionandAcceptance,nevents,colour,legend);
      if(Lumi>0){
	P.Plot1DSignificance((*Extradist1d.at(i)),true,false,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
	P.Plot1DSignificance((*Extradist1d.at(i)),false,true,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
	P.Plot1Dsigtobkg((*Extradist1d.at(i)),true,false,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
	P.Plot1Dsigtobkg((*Extradist1d.at(i)),false,true,Lumi,CrossSectionandAcceptance,nevents,colour,legend);
	P.Plot1D_DataMC_Compare((*Extradist1d.at(i)),Lumi,CrossSectionandAcceptance,nevents,colour,legend);
      }
    }
    for(unsigned int i=0; i<Extradist2d.size();i++){
      P.Plot2D((*Extradist2d.at(i)),Lumi,CrossSectionandAcceptance,nevents,colour,legend);
    }
    
    std::cout << "Writing out "<< Name << ".tex" << std::endl;
    Tables T(Name);
    T.MakeEffTable(Npassed,title,Lumi,CrossSectionandAcceptance,nevents);
    std::cout << "Plots and Tables Complete"<< std::endl;
  }
  std::cout << "H" << std::endl;
  //Check that the correct number of events are run over
  //SC.CorrectNEvents(types,nevents_noweight);
  std::cout << "I" << std::endl;
}


bool Selection::Passed(){
  for(int i=0; i<pass.size();i++){
    if(!pass.at(i)) return false;
  }
  return true;
}


double Selection::Compute(double thisdata,double thissignal, double thissignalTotal, double thisbkg, 
			  double data,double signal,double signalTotal, double bkg){

  double thisCS=(thisdata-thisbkg)*(thissignal/thissignalTotal);
  double CS=(data-bkg)*(signal/signalTotal);
  return fabs(thisCS-CS);
}

void Selection::EvaluateSystematics(Selection_Base* &selectionsys, double w){
  Selection *selsys=(Selection*)selectionsys;
  for(int j=0; j<Npassed.size();j++){
    for(int l=0; l<=Npassed.at(j).GetNbinsX();l++){
      double err=Npassed.at(j).GetBinError(l);
      if(Npassed.at(j).GetBinContent(l)!=0){
	Npassed.at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Npassed().at(j).GetBinContent(l)-Npassed.at(j).GetBinContent(l),2.0)));
      }
    }
    for(int k=0;k<Nminus1.size();k++){
      for(int l=0; l<=Nminus1.at(k).at(j).GetNbinsX();l++){
	double err=Nminus1.at(k).at(j).GetBinError(l);
	Nminus1.at(k).at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Nminus1().at(k).at(j).GetBinContent(l)-Nminus1.at(k).at(j).GetBinContent(l),2.0)));
      }
      for(int l=0; l<=Nminus0.at(k).at(j).GetNbinsX();l++){
	double err=Nminus0.at(k).at(j).GetBinError(l);
	Nminus0.at(k).at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Nminus0().at(k).at(j).GetBinContent(l)-Nminus0.at(k).at(j).GetBinContent(l),2.0)));
      }
      if(distindx.at(k)){
	for(int l=0; l<=Nminus1dist.at(k).at(j).GetNbinsX();l++){
	  double err=Nminus1dist.at(k).at(j).GetBinError(l);
	  Nminus1dist.at(k).at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Nminus1dist().at(k).at(j).GetBinContent(l)-Nminus1dist.at(k).at(j).GetBinContent(l),2.0)));
	}
	for(int l=0; l<=Accumdist.at(k).at(j).GetNbinsX();l++){
	  double err=Accumdist.at(k).at(j).GetBinError(l);
	  Accumdist.at(k).at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Accumdist().at(k).at(j).GetBinContent(l)-Accumdist.at(k).at(j).GetBinContent(l),2.0)));
	}
      }
    }
    for(int k=0;k<Extradist1d.size();k++){
      for(int l=0; l<=Extradist1d.at(k)->at(j).GetNbinsX();l++){
	double err=Extradist1d.at(k)->at(j).GetBinError(l);
	Extradist1d.at(k)->at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Extradist1d().at(k)->at(j).GetBinContent(l)-Extradist1d.at(k)->at(j).GetBinContent(l),2.0)));
      }
    }
    for(int k=0;k<Extradist2d.size();k++){
      for(int l=0; l<=Extradist2d.at(k)->at(j).GetBin(Extradist2d.at(k)->at(j).GetNbinsX(),Extradist2d.at(k)->at(j).GetNbinsY());l++){
	double err=Extradist2d.at(k)->at(j).GetBinError(l);
	Extradist2d.at(k)->at(j).SetBinError(l,sqrt(err*err+w*w*pow(selsys->Get_Extradist2d().at(k)->at(j).GetBinContent(l)-Extradist2d.at(k)->at(j).GetBinContent(l),2.0)));
      }
    }
  }
}



void Selection::ResetEvent(){
  for(int i=0; i<pass.size();i++){
    pass.at(i)=false;
  }
  for(int i=0; i<dist.size();i++){
    dist.at(i).clear();
  }
}