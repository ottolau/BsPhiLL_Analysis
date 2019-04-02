#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TRandom.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"

typedef struct {
  Float_t       m_ee;
  Float_t       m_KK;
  Float_t       m_KKee;
  Float_t       mvaValue;
  Float_t	eledR;
  Float_t	kaondR;
  Float_t	jpsiPhidR;
} mvaTree_t;

void BsPhiLLTupleTree::Loop(TString outputfile, Int_t maxevents, Float_t mvaCut, bool oppCharge)
{
//   In a ROOT session, you can do:
//      root> .L BsPhiLLTupleTree.C
//      root> BsPhiLLTupleTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
    fChain->SetBranchStatus("*",0);  // disable all branches
    fChain->SetBranchStatus("nEle",1);  // activate branchname
    fChain->SetBranchStatus("eleCharge_lead",1);
    fChain->SetBranchStatus("eleD0_lead",1);
    fChain->SetBranchStatus("eleDz_lead",1);
    fChain->SetBranchStatus("eleD0Error_lead",1);
    fChain->SetBranchStatus("eleDzError_lead",1);
    fChain->SetBranchStatus("elePt_lead",1);
    fChain->SetBranchStatus("eleEta_lead",1);
    fChain->SetBranchStatus("elePhi_lead",1);
    fChain->SetBranchStatus("eleConvVeto_lead",1);
    fChain->SetBranchStatus("eleCharge_sublead",1);
    fChain->SetBranchStatus("eleD0_sublead",1);
    fChain->SetBranchStatus("eleDz_sublead",1);
    fChain->SetBranchStatus("eleD0Error_sublead",1);
    fChain->SetBranchStatus("eleDzError_sublead",1);
    fChain->SetBranchStatus("elePt_sublead",1);
    fChain->SetBranchStatus("eleEta_sublead",1);
    fChain->SetBranchStatus("elePhi_sublead",1);
    fChain->SetBranchStatus("eleConvVeto_sublead",1);
    fChain->SetBranchStatus("eleSvChi2",1);
    fChain->SetBranchStatus("eleSvNDOF",1);
    fChain->SetBranchStatus("eleSvProb",1);
    fChain->SetBranchStatus("eleSvCtxy",1);
    fChain->SetBranchStatus("eleSvCosAngle",1);
    fChain->SetBranchStatus("eleSvLxy",1);
    fChain->SetBranchStatus("eleSvLxyError",1);
    fChain->SetBranchStatus("kaonEECharge_lead",1);
    fChain->SetBranchStatus("kaonEED0_lead",1);
    fChain->SetBranchStatus("kaonEEDz_lead",1);
    fChain->SetBranchStatus("kaonEED0Error_lead",1);
    fChain->SetBranchStatus("kaonEEDzError_lead",1);
    fChain->SetBranchStatus("kaonEEPt_lead",1);
    fChain->SetBranchStatus("kaonEEEta_lead",1);
    fChain->SetBranchStatus("kaonEEPhi_lead",1);
    fChain->SetBranchStatus("kaonEETrkNormChi2_lead",1);
    fChain->SetBranchStatus("kaonEECharge_sublead",1);
    fChain->SetBranchStatus("kaonEED0_sublead",1);
    fChain->SetBranchStatus("kaonEEDz_sublead",1);
    fChain->SetBranchStatus("kaonEED0Error_sublead",1);
    fChain->SetBranchStatus("kaonEEDzError_sublead",1);
    fChain->SetBranchStatus("kaonEEPt_sublead",1);
    fChain->SetBranchStatus("kaonEEEta_sublead",1);
    fChain->SetBranchStatus("kaonEEPhi_sublead",1);
    fChain->SetBranchStatus("kaonEETrkNormChi2_sublead",1);
    fChain->SetBranchStatus("bsEEdRele",1);
    fChain->SetBranchStatus("bsEEdRkaon",1);
    fChain->SetBranchStatus("bsEEdRJpsiPhi",1);
    fChain->SetBranchStatus("bsEEJpsiMass",1);
    fChain->SetBranchStatus("bsEEPhiMass",1);
    fChain->SetBranchStatus("bsEEBsMass",1);

// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile file_out(outputfile, "recreate");

   // initialization
   Float_t       elePtLead = 0;
   Float_t       eleD0Lead = 0;
   Float_t       eleDzLead = 0;
   Float_t       elePtSublead = 0;
   Float_t       eleD0Sublead =  0;
   Float_t       eleDzSublead = 0;
   Float_t       kaonPtLead = 0;
   Float_t       kaonD0Lead = 0;
   Float_t       kaonDzLead = 0;
   Float_t       kaonNormChi2Lead = 0;
   Float_t       kaonPtSublead = 0;
   Float_t       kaonD0Sublead = 0;
   Float_t       kaonDzSublead = 0;
   Float_t       kaonNormChi2Sublead = 0;
   Float_t       eledR = 0;
   Float_t       kaondR = 0;
   Float_t       jpsiPhidR = 0;
   Float_t       svProb = 0;
   Float_t       svCosine = 0;
   Float_t       svLxySig = 0;
   Float_t       jpsiPt = 0;
   Float_t       phiPt = 0;
   Float_t       bsPt = 0;
   Float_t       phiMassFrac = 0;

   TMVA::Tools::Instance();
   TMVA::PyMethodBase::PyInitialize();

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   reader->AddVariable("elePtLead", &elePtLead);
   reader->AddVariable("elePtSublead", &elePtSublead);
   reader->AddVariable("kaonPtLead", &kaonPtLead);
   reader->AddVariable("kaonPtSublead", &kaonPtSublead);
   //reader->AddVariable("jpsiPt", &jpsiPt);
   //reader->AddVariable("phiPt", &phiPt);
   reader->AddVariable("bsPt", &bsPt);
   //reader->AddVariable("eledR", &eledR);
   //reader->AddVariable("kaondR", &kaondR);
   //reader->AddVariable("jpsiPhidR", &jpsiPhidR);
   reader->AddVariable("svProb", &svProb);
   reader->AddVariable("svCosine", &svCosine);
   reader->AddVariable("svLxySig", &svLxySig);
   reader->AddVariable("eleD0Lead", &eleD0Lead);
   reader->AddVariable("eleD0Sublead", &eleD0Sublead);
   reader->AddVariable("eleDzLead", &eleDzLead);
   reader->AddVariable("eleDzSublead", &eleDzSublead);
   reader->AddVariable("kaonD0Lead", &kaonD0Lead);
   reader->AddVariable("kaonD0Sublead", &kaonD0Sublead);
   reader->AddVariable("kaonDzLead", &kaonDzLead);
   reader->AddVariable("kaonDzSublead", &kaonDzSublead);
   reader->AddVariable("kaonNormChi2Lead", &kaonNormChi2Lead);
   reader->AddVariable("kaonNormChi2Sublead", &kaonNormChi2Sublead);
   //reader->AddVariable("phiMassFrac", &phiMassFrac);

   TString TMVAmethod = "PyKeras method";

   reader->BookMVA("PyKeras method", "/uscms/home/klau/nobackup/BtoKll/BParking2018/BsPhiLL_Integrated/BsPhiLL_Analysis/dataset/weights/TMVAClassification_BsPhiJpsiEE_noPhiM_pTcuts_PyKeras.weights.xml");

   // Initialization
   TTree *tree = new TTree("mvaTree", "mva tree");
   mvaTree_t mvaTree;
   tree->Branch("m_ee",      &mvaTree.m_ee,      "m_ee/F");
   tree->Branch("m_KK",      &mvaTree.m_KK,      "m_KK/F");
   tree->Branch("m_KKee",    &mvaTree.m_KKee,    "m_KKee/F");
   tree->Branch("mvaValue",  &mvaTree.mvaValue,  "mvaValue/F");
   tree->Branch("eledR",     &mvaTree.eledR,     "eledR/F");
   tree->Branch("kaondR",    &mvaTree.kaondR,    "kaondR/F");
   tree->Branch("jpsiPhidR", &mvaTree.jpsiPhidR, "jpsiPhidR/F");

   Int_t nAccept = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (nAccept >= maxevents && maxevents != -1) break;    
      if ((jentry % 10000) == 0) cout<<"Processing entry: "<<jentry<<", accepted events: "<<nAccept<<endl;

      /***************************
       *
       * Electron Part
       *
       * ************************/

      double bestMVA = -1.0;
      double bestMVAarg = -1;

      for (int i=0; i<nEle; ++i) {

	 bool preCut = true;
	
	 if (oppCharge) {
	   if (eleCharge_lead->at(i) * eleCharge_sublead->at(i) > 0) continue;
	   if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) > 0) continue;
	 } else {
	   if (eleCharge_lead->at(i) * eleCharge_sublead->at(i) < 0) continue;
	   if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) < 0) continue;
	 }
	 if (elePt_lead->at(i) < 0.4 || elePt_sublead->at(i) < 0.4 || kaonEEPt_lead->at(i) < 0.4 || kaonEEPt_sublead->at(i) < 0.4) continue;
         if (!eleConvVeto_lead->at(i) || !eleConvVeto_sublead->at(i)) continue;
	 if (bsEEJpsiMass->at(i) < jpsilow || bsEEJpsiMass->at(i) > jpsiup) continue;
	 if (bsEEPhiMass->at(i) < philow || bsEEPhiMass->at(i) > phiup) continue;

	 elePtLead = elePt_lead->at(i)/bsEEBsMass->at(i);
	 eleD0Lead = eleD0_lead->at(i);
	 eleDzLead = eleDz_lead->at(i);
	 elePtSublead = elePt_sublead->at(i)/bsEEBsMass->at(i);
	 eleD0Sublead = eleD0_sublead->at(i);
	 eleDzSublead = eleDz_sublead->at(i);
	 kaonPtLead = kaonEEPt_lead->at(i)/bsEEBsMass->at(i);
	 kaonD0Lead = kaonEED0_lead->at(i);
	 kaonDzLead = kaonEEDz_lead->at(i);
	 kaonNormChi2Lead = kaonEETrkNormChi2_lead->at(i);
	 kaonPtSublead = kaonEEPt_sublead->at(i)/bsEEBsMass->at(i);
	 kaonD0Sublead = kaonEED0_sublead->at(i);
	 kaonDzSublead = kaonEEDz_sublead->at(i);
	 kaonNormChi2Sublead = kaonEETrkNormChi2_sublead->at(i);
	 eledR = bsEEdRele->at(i);
	 kaondR = bsEEdRkaon->at(i);
	 jpsiPhidR = bsEEdRJpsiPhi->at(i);
	 jpsiPt = TMath::Sqrt(pow(elePt_lead->at(i),2) + pow(elePt_sublead->at(i),2) + 2.0*elePt_lead->at(i)*elePt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),elePhi_sublead->at(i))))/bsEEBsMass->at(i);
	 phiPt = TMath::Sqrt(pow(kaonEEPt_lead->at(i),2) + pow(kaonEEPt_sublead->at(i),2) + 2.0*kaonEEPt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonEEPhi_lead->at(i),kaonEEPhi_sublead->at(i))))/bsEEBsMass->at(i);
	 bsPt = TMath::Sqrt(pow(jpsiPt*bsEEBsMass->at(i),2) + pow(phiPt*bsEEBsMass->at(i),2) 
	  + 2.0*elePt_lead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_lead->at(i))) 
	  + 2.0*elePt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_sublead->at(i))) 
	  + 2.0*elePt_sublead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_lead->at(i))) 
	  + 2.0*elePt_sublead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_sublead->at(i))) 
	  )/bsEEBsMass->at(i);
	 svProb = eleSvProb->at(i);
	 svCosine = eleSvCosAngle->at(i);
	 svLxySig = eleSvLxy->at(i)/eleSvLxyError->at(i);
	 phiMassFrac = bsEEPhiMass->at(i)/phiM;

	 double mvaValue = reader->EvaluateMVA(TMVAmethod);

	 if (mvaValue > mvaCut && mvaValue > bestMVA) {
	   bestMVA = mvaValue;
	   bestMVAarg = i;
	 }
      }

      if (bestMVA > 0.0) {
	 mvaTree.m_ee = bsEEJpsiMass->at(bestMVAarg);
	 mvaTree.m_KK = bsEEPhiMass->at(bestMVAarg);
	 mvaTree.m_KKee = bsEEBsMass->at(bestMVAarg);
	 mvaTree.mvaValue = bestMVA;
	 mvaTree.eledR = bsEEdRele->at(bestMVAarg);
	 mvaTree.kaondR = bsEEdRkaon->at(bestMVAarg);
	 mvaTree.jpsiPhidR = bsEEdRJpsiPhi->at(bestMVAarg);

	 tree->Fill();
	 nAccept++;
      }
   }

   file_out.cd();
   file_out.Write();
   file_out.Close();


}
