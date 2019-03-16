#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Types.h"
#include "TMVA/CrossValidation.h"


void TMVA_DNN_Classification_BsToJpsiPhiuu_CrossValidation() {
   gROOT->SetBatch(kTRUE);
   gROOT->SetStyle("Plain");
   gStyle->SetGridStyle(3);

   TMVA::Tools::Instance();

   auto outputFile = TFile::Open("DNN_ClassificationOutput_nontrg_21params_newTree_cv.root", "RECREATE");
   TString inputFileName_bkg = "mvaTrainingSet_muon_nontrigger_bkg_Run2018ABD.root";
   TString inputFileName_sig = "mvaTrainingSet_muon_nontrigger_sig.root";

   //bool inc4vectors = true;

//   TMVA::Factory factory("TMVA_DNN_Classification_trgmu_13params", outputFile,
//	                       "!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 

//   TMVA::Factory factory("TMVA_DNN_Classification_nontrg_21params_newTree", outputFile,
//	                       "!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 


   TMVA::DataLoader * loader = new TMVA::DataLoader("dataset");

   loader->AddVariable("muPt1st");
   loader->AddVariable("muPt2nd");
   loader->AddVariable("kaonPt1st");
   loader->AddVariable("kaonPt2nd");
   loader->AddVariable("jpsiPt");
   loader->AddVariable("phiPt");
   loader->AddVariable("bsPt");
   loader->AddVariable("mudR");
   loader->AddVariable("kaondR");
   loader->AddVariable("jpsiPhidR");
   loader->AddVariable("svCtxy");
   loader->AddVariable("svChi2");
   loader->AddVariable("svCosine");
   loader->AddVariable("muIDCutBased1st");
   loader->AddVariable("muIDCutBased2nd");
   loader->AddVariable("muDxySig1st"); 
   loader->AddVariable("muDxySig2nd"); 
   loader->AddVariable("kaonDxySig1st"); 
   loader->AddVariable("kaonDxySig2nd"); 
   loader->AddVariable("kaonChi2NDF1st"); 
   loader->AddVariable("kaonChi2NDF2nd"); 


   auto inputFile_bkg = TFile::Open( inputFileName_bkg );
   auto inputFile_sig = TFile::Open( inputFileName_sig );

   TTree *signalTree     = (TTree*)inputFile_sig->Get("signal");
   TTree *backgroundTree = (TTree*)inputFile_bkg->Get("background");

   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
      
   loader->AddSignalTree    ( signalTree,     signalWeight     );
   loader->AddBackgroundTree( backgroundTree, backgroundWeight );

   TCut mycuts = ""; 
   TCut mycutb = "";

   loader->PrepareTrainingAndTestTree( mycuts, mycutb,
					     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );


   UInt_t numFolds = 5;
   TString analysisType = "Classification";
   TString splitType = "Random";
   TString splitExpr = "";

   TString cvOptions = Form("!V"
                            ":!Silent"
                            ":ModelPersistence"
                            ":AnalysisType=%s"
                            //":SplitType=%s"
                            ":NumFolds=%i"
                            ":SplitExpr=%s",
                            analysisType.Data(), /*splitType.Data(),*/ numFolds,
                            splitExpr.Data());

   TMVA::CrossValidation cv{"TMVA_DNN_Classification_nontrg_21params_CrossValidation", loader, outputFile, cvOptions};


   //Boosted Decision Trees
//   cv.BookMethod(TMVA::Types::kBDT, "BDT",
//                  "!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


   TString inputLayoutString;
   TString batchLayoutString;
   inputLayoutString = "InputLayout=1|1|21"; 
   batchLayoutString = "BatchLayout=1|128|21";
   //TString layoutString ("Layout=DENSE|32|TANH,DENSE|32|TANH,DENSE|32|TANH,DENSE|32|TANH,DENSE|1|LINEAR");
   TString layoutString ("Layout=DENSE|64|TANH,DENSE|64|TANH,DENSE|64|TANH,DENSE|64|TANH,DENSE|1|LINEAR");
   //TString layoutString ("Layout=DENSE|128|TANH,DENSE|128|TANH,DENSE|128|TANH,DENSE|128|TANH,DENSE|1|LINEAR");
   //TString layoutString ("Layout=DENSE|8|TANH,DENSE|8|TANH,DENSE|8|TANH,DENSE|8|TANH,DENSE|1|LINEAR");
																					    
   //Training strategy
   TString training1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
		     "ConvergenceSteps=20,BatchSize=128,TestRepetitions=1,"
		     "MaxEpochs=300,WeightDecay=1e-4,Regularization=L2,"
		     "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.");

   TString training2("LearningRate=1e-4,Momentum=0.9,Repetitions=1,"
		     "ConvergenceSteps=20,BatchSize=128,TestRepetitions=1,"
		     "MaxEpochs=50,WeightDecay=1e-4,Regularization=L2,"
		     "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.");

   TString trainingStrategyString ("TrainingStrategy=");
   trainingStrategyString += training1 + "|" + training2; //+ "|" + training3;

   //Options                                                                                                                                                                
   TString dnnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
		       "WeightInitialization=XAVIERUNIFORM");
   dnnOptions.Append (":"); dnnOptions.Append (inputLayoutString);
   dnnOptions.Append (":"); dnnOptions.Append (batchLayoutString);
   dnnOptions.Append (":"); dnnOptions.Append (layoutString);
   dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

   dnnOptions += ":Architecture=Standard";
   cv.BookMethod(TMVA::Types::kDL, "DL_DENSE", dnnOptions);

   // --------------------------------------------------------------------------

   //
   // Train, test and evaluate the booked methods.
   // Evaluates the booked methods once for each fold and aggregates the result
   // in the specified output file.
   //
   cv.Evaluate();

   // --------------------------------------------------------------------------

   //
   // Process some output programatically, printing the ROC score for each
   // booked method.
   //
   size_t iMethod = 0;
   for (auto && result : cv.GetResults()) {
      std::cout << "Summary for method " << cv.GetMethods()[iMethod++].GetValue<TString>("MethodName")
                << std::endl;
      for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
         std::cout << "\tFold " << iFold << ": "
                   << "ROC int: " << result.GetROCValues()[iFold]
                   << ", "
                   << "BkgEff@SigEff=0.3: " << result.GetEff30Values()[iFold]
                   << std::endl;
      }
   }

   // --------------------------------------------------------------------------

   //
   // Save the output
   //
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVACrossValidation is done!" << std::endl;

   // --------------------------------------------------------------------------



}


