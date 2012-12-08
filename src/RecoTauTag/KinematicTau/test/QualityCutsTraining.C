// @(#)root/tmva $Id: TMVAClassification.C 31458 2009-11-30 13:58:20Z stelzer $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l root -l QualityCutsTraining.C\(\"BDT,Likelihood\",\"+\"\)           *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/


#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void QualityCutsTraining(TString myMethodList = "", string type = "+"){


  TMVA::Tools::Instance();
  std::map<std::string,int> Use;


   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 1;
   Use["CutsGA"]          = 1;
   Use["CutsSA"]          = 1;
   // ---
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 1; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 1;
   Use["LikelihoodMIX"]   = 1;
   // ---
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 1;
   Use["PDERSPCA"]        = 1;
   Use["PDERSkNN"]        = 1; // depreciated until further notice
   Use["PDEFoam"]         = 1;
   // --
   Use["KNN"]             = 1;
   // ---
   Use["HMatrix"]         = 1;
   Use["Fisher"]          = 1;
   Use["FisherG"]         = 1;
   Use["BoostedFisher"]   = 1;
   Use["LD"]              = 1;
   // ---
   Use["FDA_GA"]          = 1;
   Use["FDA_SA"]          = 1;
   Use["FDA_MC"]          = 1;
   Use["FDA_MT"]          = 1;
   Use["FDA_GAMT"]        = 1;
   Use["FDA_MCMT"]        = 1;
   // ---
   Use["MLP"]             = 1; // this is the recommended ANN
   Use["MLPBFGS"]         = 1; // recommended ANN with optional training method
   Use["CFMlpANN"]        = 1; // *** missing
   Use["TMlpANN"]         = 1; 
   // ---
   Use["SVM"]             = 1;
   // ---
   Use["BDT"]             = 1;
   Use["BDTD"]            = 1;
   Use["BDTG"]            = 1;
   Use["BDTB"]            = 1;
   // ---
   Use["RuleFit"]         = 1;
   // ---
   Use["Plugin"]          = 0;
   // ---------------------------------------------------------------


   std::cout << std::endl;
   std::cout << "==> Start QualityCutsTraining" << std::endl;
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // Create a new root output file.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "QualityCutsTraining", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );





   //  TFile * dy = TFile::Open("/home/home2/institut_3b/cherepanov/work/KinFitDevelop/CMSSW_5_2_5/src/DTBasedCutsOptimization/output/dy/v2/dy.root");
   TFile * dy = TFile::Open("/user/cherepanov/SkimForBDT/DY_v8/dy.root");
   //   TFile * wl = TFile::Open("/home/home2/institut_3b/cherepanov/work/KinFitDevelop/CMSSW_5_2_5/src/DTBasedCutsOptimization/output/wl/v2/wl.root");
   //   TFile * qcd = TFile::Open("/home/home2/institut_3b/cherepanov/work/KinFitDevelop/CMSSW_5_2_5/src/DTBasedCutsOptimization/output/QCD/QCD.root");
   TFile * qcd = TFile::Open("/user/cherepanov/SkimForBDT/QCD/qcd.root");

	float fracPlus;
	float a1MassPlus;
	float ProbPlus2;
	float iterPlus;
	float PVSVPlus;
	float PVPVPlus;



	float fracMins;
	float a1MassMins;
	float ProbMins2;
	float iterMins;
	float PVSVMins;
	float PVPVMins;



	float fracZero;
	float a1MassZero;
	float ProbZero2;
	float iterZero;
	float PVSVZero;
	float PVPVZero;


  //TMVA_dy.root
	if(type == "+"){
	  factory->AddVariable( "fracPlus", 'F' );
	  factory->AddVariable( "a1MassPlus", 'F' );
	  factory->AddVariable( "ProbPlus3", 'F' );
	  factory->AddVariable( "iterPlus", 'F' );
	  factory->AddVariable( "PVSVPlus", 'F' );
	}else  if(type == "-"){
	  factory->AddVariable( "fracMins", 'F' );
	  factory->AddVariable( "a1MassMins", 'F' );
	  factory->AddVariable( "ProbMins3", 'F' );
	  factory->AddVariable( "iterMins", 'F' );
	  factory->AddVariable( "PVSVMins", 'F' );
	}else  if(type=="0"){
	  factory->AddVariable( "fracZero", 'F' );
	  factory->AddVariable( "a1MassZero", 'F' );
	  factory->AddVariable( "ProbZero3", 'F' );
	  factory->AddVariable( "PVSVZero", 'F' );
	  factory->AddVariable( "PVPVZero", 'F' );
	}
	

      TTree *signal      = (TTree*)dy->Get("TrainTree");
      TTree *background1 = (TTree*)qcd->Get("TrainTree");


      // global event weights per tree (see below for setting event-wise weights)
      Double_t signalWeight     = 1.0;
      Double_t backgroundWeight1 = 1.0;


      factory->AddSignalTree( signal,     signalWeight     );
      factory->AddBackgroundTree( background1, backgroundWeight1 );


  //  if ( vtxSignPVRotSV < 2. )return false; // Sig. of secondary vertex
  //  if ( vtxSignPVRotPVRed > 2. )return false; //vertex sig. between modified and initial primary vertex
  //WARNING!!!
  //from now one we assume a tau decay into three pions and neutrino
  //other channels need their own discriminators
  //!!!

  //  if(chargedDaughters.size()!=1 || neutralDaughters.size()!=1) return false; // number of decay products
  // if(KFTau.PFTauRef()->signalPFChargedHadrCands().size() > 3 ) return false; //tracks in signal cone of initial pftau candidate
  //  if(a1Mass < 0.8) return false; //refitPFTau equals refitted a1 in 3-prong case
  //  if(energyTFraction == -1.) return false; //energy fraction
  // if(energyTFraction < 0 || energyTFraction > 1.) return false;  //energy fraction


      if(type == "+"){
	TCut mycuts = "fracPlus < 1 && fracPlus >0 && a1MassPlus > 0.8 && PVSVPlus > 2 && ProbPlus3>0.03";
      }else if(type == "-"){
	//  TCut mycuts = "";
		TCut mycuts = "fracMins < 0.9 && fracMins >0.6 && a1MassMins > 0.85 &&  PVSVMins > 2 && ProbMins3>0.03";
      }else if(type == "0"){
	TCut mycuts = "fracZero < 0.9 && fracZero >0.6 && a1MassZero > 0.85  && PVSVZero > 2 && PVPVZero <2 && ProbZero3>0.03";
      }
      TCut mycutb = ""; 
      
      factory->PrepareTrainingAndTestTree( mycuts, mycutb,
					   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

      if (Use["Likelihood"])
	factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", 
			     "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" ); 
      if (Use["BDT"])  // Adaptive Boost
	factory->BookMethod( TMVA::Types::kBDT, "BDT", 
			     //                          "!H:!V:NTrees=1000:nEventsMin=x1:MaxDepth=1000:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=-1:PruneMethod=CostComplexity" );
                             "!H:!V:NTrees=100:nEventsMin=2000:MaxDepth=500:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=-1:PruneMethod=NoPruning" );
 



   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> |TrainPOlarization is done!" << std::endl;      

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
 
