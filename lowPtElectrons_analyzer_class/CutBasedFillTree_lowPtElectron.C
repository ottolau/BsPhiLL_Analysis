#define LowPtElectrons_BsPhiLLTupleTree_cxx
#include "LowPtElectrons_BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"

typedef struct {
  Float_t       m_ee;
  Float_t       m_KK;
  Float_t       m_KKee;
  Float_t	elePtLead;
  Float_t	elePtSublead;
  Float_t	kaonPtLead;
  Float_t	kaonPtSublead;
  Float_t       bsPt;
  Float_t       LxySig;
  Float_t       SvProb;
  Float_t       SvCosine;
  Float_t	mvaBWP_lead;
  Float_t	mvaUnBWP_lead;
  Float_t	mvaBWP_sublead;
  Float_t	mvaUnBWP_sublead;
  Float_t	eledR;
  Float_t	kaondR;
  Float_t	jpsiPhidR;
} cutTree_t;

void LowPtElectrons_BsPhiLLTupleTree::Loop(TString outputfile, Int_t maxevents, Float_t mvaCut, bool oppCharge)
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
    fChain->SetBranchStatus("nLowPt",1);  // activate branchname
    fChain->SetBranchStatus("lowPtCharge_lead",1);
    fChain->SetBranchStatus("lowPtD0_lead",1);
    fChain->SetBranchStatus("lowPtDz_lead",1);
    fChain->SetBranchStatus("lowPtD0Error_lead",1);
    fChain->SetBranchStatus("lowPtDzError_lead",1);
    fChain->SetBranchStatus("lowPtPt_lead",1);
    fChain->SetBranchStatus("lowPtEta_lead",1);
    fChain->SetBranchStatus("lowPtPhi_lead",1);
    fChain->SetBranchStatus("lowPtMVABWP_lead",1);
    fChain->SetBranchStatus("lowPtMVAUnBWP_lead",1);
    fChain->SetBranchStatus("lowPtCharge_sublead",1);
    fChain->SetBranchStatus("lowPtD0_sublead",1);
    fChain->SetBranchStatus("lowPtDz_sublead",1);
    fChain->SetBranchStatus("lowPtD0Error_sublead",1);
    fChain->SetBranchStatus("lowPtDzError_sublead",1);
    fChain->SetBranchStatus("lowPtPt_sublead",1);
    fChain->SetBranchStatus("lowPtEta_sublead",1);
    fChain->SetBranchStatus("lowPtPhi_sublead",1);
    fChain->SetBranchStatus("lowPtMVABWP_sublead",1);
    fChain->SetBranchStatus("lowPtMVAUnBWP_sublead",1);
    fChain->SetBranchStatus("lowPtSvChi2",1);
    fChain->SetBranchStatus("lowPtSvNDOF",1);
    fChain->SetBranchStatus("lowPtSvProb",1);
    fChain->SetBranchStatus("lowPtSvCtxy",1);
    fChain->SetBranchStatus("lowPtSvCosAngle",1);
    fChain->SetBranchStatus("lowPtSvLxy",1);
    fChain->SetBranchStatus("lowPtSvLxyError",1);
    fChain->SetBranchStatus("kaonLowPtCharge_lead",1);
    fChain->SetBranchStatus("kaonLowPtD0_lead",1);
    fChain->SetBranchStatus("kaonLowPtDz_lead",1);
    fChain->SetBranchStatus("kaonLowPtD0Error_lead",1);
    fChain->SetBranchStatus("kaonLowPtDzError_lead",1);
    fChain->SetBranchStatus("kaonLowPtPt_lead",1);
    fChain->SetBranchStatus("kaonLowPtEta_lead",1);
    fChain->SetBranchStatus("kaonLowPtPhi_lead",1);
    fChain->SetBranchStatus("kaonLowPtTrkNormChi2_lead",1);
    fChain->SetBranchStatus("kaonLowPtCharge_sublead",1);
    fChain->SetBranchStatus("kaonLowPtD0_sublead",1);
    fChain->SetBranchStatus("kaonLowPtDz_sublead",1);
    fChain->SetBranchStatus("kaonLowPtD0Error_sublead",1);
    fChain->SetBranchStatus("kaonLowPtDzError_sublead",1);
    fChain->SetBranchStatus("kaonLowPtPt_sublead",1);
    fChain->SetBranchStatus("kaonLowPtEta_sublead",1);
    fChain->SetBranchStatus("kaonLowPtPhi_sublead",1);
    fChain->SetBranchStatus("kaonLowPtTrkNormChi2_sublead",1);
    fChain->SetBranchStatus("bsLowPtdRele",1);
    fChain->SetBranchStatus("bsLowPtdRkaon",1);
    fChain->SetBranchStatus("bsLowPtdRJpsiPhi",1);
    fChain->SetBranchStatus("bsLowPtJpsiMass",1);
    fChain->SetBranchStatus("bsLowPtPhiMass",1);
    fChain->SetBranchStatus("bsLowPtBsMass",1);


// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile file_out(outputfile, "recreate");

   // Initialization
   TTree *tree = new TTree("cutTree", "cut-based tree");
   cutTree_t cutTree;
   tree->Branch("m_ee",      &cutTree.m_ee,      "m_ee/F");
   tree->Branch("m_KK",      &cutTree.m_KK,      "m_KK/F");
   tree->Branch("m_KKee",    &cutTree.m_KKee,    "m_KKee/F");
   tree->Branch("elePtLead", &cutTree.elePtLead, "elePtLead/F");
   tree->Branch("elePtSublead", &cutTree.elePtSublead, "elePtSublead/F");
   tree->Branch("kaonPtLead",&cutTree.kaonPtLead,"kaonPtLead/F");
   tree->Branch("kaonPtSublead", &cutTree.kaonPtSublead, "kaonPtSublead/F");
   tree->Branch("bsPt",      &cutTree.bsPt,      "bsPt/F");
   tree->Branch("LxySig",    &cutTree.LxySig,    "LxySig/F");
   tree->Branch("SvProb",    &cutTree.SvProb,    "SvProb/F");
   tree->Branch("SvCosine",  &cutTree.SvCosine,  "SvCosine/F");
   tree->Branch("mvaBWP_lead",    &cutTree.mvaBWP_lead,    "mvaBWP_lead/F");
   tree->Branch("mvaUnBWP_lead",  &cutTree.mvaUnBWP_lead,  "mvaUnBWP_lead/F");
   tree->Branch("mvaBWP_sublead", &cutTree.mvaBWP_sublead, "mvaBWP_sublead/F");
   tree->Branch("mvaUnBWP_sublead",  &cutTree.mvaUnBWP_sublead,  "mvaUnBWP_sublead/F");
   tree->Branch("eledR",     &cutTree.eledR,     "eledR/F");
   tree->Branch("kaondR",    &cutTree.kaondR,    "kaondR/F");
   tree->Branch("jpsiPhidR", &cutTree.jpsiPhidR, "jpsiPhidR/F");


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if ((jentry%100000) == 0) std::cout<<"Procssing entry: "<<jentry<<std::endl;


      /***************************
       *
       * Electron Part
       *
       * ************************/


      for (int i=0; i<nLowPt; ++i) {

	 bool preCut = true;

	 if (lowPtCharge_lead->at(i) * lowPtCharge_sublead->at(i) < 0) continue;
	 if (kaonLowPtCharge_lead->at(i) * kaonLowPtCharge_sublead->at(i) > 0) continue; // same charge
	 //if (lowPtPt_lead->at(i) < 0.4 || lowPtPt_sublead->at(i) < 0.4 || kaonLowPtPt_lead->at(i) < 0.4 || kaonLowPtPt_sublead->at(i) < 0.4) continue;

	 bool cutBasedSelected = true;
	 //if (lowPtPt_lead->at(i) < 2 || lowPtPt_sublead->at(i) < 2) cutBasedSelected = false;
	 //if (kaonLowPtPt_lead->at(i) < 0.8 || kaonLowPtPt_sublead->at(i) < 0.8) cutBasedSelected = false;
	 if (lowPtSvProb->at(i) < 0.01) cutBasedSelected = false;
	 //if (lowPtSvLxy->at(i)/lowPtSvLxyError->at(i) < 3) cutBasedSelected = false;
	 if (lowPtSvCosAngle->at(i) < 0.9) cutBasedSelected = false;
	 //if (bsLowPtdRele->at(i) > 1.5) cutBasedSelected = false;
	 //if (bsLowPtdRkaon->at(i) > 0.5) cutBasedSelected = false;
	 //if (bsLowPtdRJpsiPhi->at(i) > 1.5) cutBasedSelected = false;
	 //if (lowPtMVAUnBWP_lead->at(i) < 6.0) cutBasedSelected = false;
	 if (lowPtMVAUnBWP_lead->at(i) < 3.05) cutBasedSelected = false;
	 if (lowPtMVAUnBWP_sublead->at(i) < 3.05) cutBasedSelected = false;

	 double jpsiPt = TMath::Sqrt(pow(lowPtPt_lead->at(i),2) + pow(lowPtPt_sublead->at(i),2) + 2.0*lowPtPt_lead->at(i)*lowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),lowPtPhi_sublead->at(i))))/bsLowPtBsMass->at(i);
	 double phiPt = TMath::Sqrt(pow(kaonLowPtPt_lead->at(i),2) + pow(kaonLowPtPt_sublead->at(i),2) + 2.0*kaonLowPtPt_lead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonLowPtPhi_lead->at(i),kaonLowPtPhi_sublead->at(i))))/bsLowPtBsMass->at(i);
	 double bsPt = TMath::Sqrt(pow(jpsiPt*bsLowPtBsMass->at(i),2) + pow(phiPt*bsLowPtBsMass->at(i),2) 
		     + 2.0*lowPtPt_lead->at(i)*kaonLowPtPt_lead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),kaonLowPtPhi_lead->at(i))) 
		     + 2.0*lowPtPt_lead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),kaonLowPtPhi_sublead->at(i))) 
		     + 2.0*lowPtPt_sublead->at(i)*kaonLowPtPt_lead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_sublead->at(i),kaonLowPtPhi_lead->at(i))) 
		     + 2.0*lowPtPt_sublead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_sublead->at(i),kaonLowPtPhi_sublead->at(i))) 
		     )/bsLowPtBsMass->at(i);

	 //if (bsPt < 3) cutBasedSelected = false;

	 if (cutBasedSelected) {

	    cutTree.m_ee = bsLowPtJpsiMass->at(i);
	    cutTree.m_KK = bsLowPtPhiMass->at(i);
	    cutTree.m_KKee = bsLowPtBsMass->at(i);
	    cutTree.elePtLead = lowPtPt_lead->at(i);
	    cutTree.elePtSublead = lowPtPt_sublead->at(i);
	    cutTree.kaonPtLead = kaonLowPtPt_lead->at(i);
	    cutTree.kaonPtSublead = kaonLowPtPt_sublead->at(i);
	    cutTree.bsPt = bsPt*bsLowPtBsMass->at(i);
	    cutTree.LxySig = lowPtSvLxy->at(i)/lowPtSvLxyError->at(i);
	    cutTree.SvProb = lowPtSvProb->at(i);
	    cutTree.SvCosine = lowPtSvCosAngle->at(i);
	    cutTree.mvaBWP_lead = lowPtMVABWP_lead->at(i);
	    cutTree.mvaUnBWP_lead = lowPtMVAUnBWP_lead->at(i);
	    cutTree.mvaBWP_sublead = lowPtMVABWP_sublead->at(i);
	    cutTree.mvaUnBWP_sublead = lowPtMVAUnBWP_sublead->at(i);
	    cutTree.eledR = bsLowPtdRele->at(i);
	    cutTree.kaondR = bsLowPtdRkaon->at(i);
	    cutTree.jpsiPhidR = bsLowPtdRJpsiPhi->at(i);

	    tree->Fill();


	 }

      }


   }

   file_out.cd();
   file_out.Write();
   file_out.Close();

}
