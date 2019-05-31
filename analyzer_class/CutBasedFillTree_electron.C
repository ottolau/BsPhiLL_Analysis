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
  Float_t	elePtLead;
  Float_t	elePtSublead;
  Float_t	kaonPtLead;
  Float_t	kaonPtSublead;
  Float_t       bsPt;
  Float_t       LxySig;
  Float_t       SvProb;
  Float_t       SvCosine;
  Float_t	eledR;
  Float_t	kaondR;
  Float_t	jpsiPhidR;
  Float_t       eleID_lead;
  Float_t	eleID_sublead;
} cutTree_t;

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
    fChain->SetBranchStatus("elePt_lead",1);
    fChain->SetBranchStatus("eleEta_lead",1);
    fChain->SetBranchStatus("elePhi_lead",1);
    fChain->SetBranchStatus("eleIDbit_lead",1);
    fChain->SetBranchStatus("eleConvVeto_lead",1);
    fChain->SetBranchStatus("eleCharge_sublead",1);
    fChain->SetBranchStatus("elePt_sublead",1);
    fChain->SetBranchStatus("eleEta_sublead",1);
    fChain->SetBranchStatus("elePhi_sublead",1);
    fChain->SetBranchStatus("eleIDbit_sublead",1);
    fChain->SetBranchStatus("eleConvVeto_sublead",1);
    fChain->SetBranchStatus("eleSvProb",1);
    fChain->SetBranchStatus("eleSvCtxy",1);
    fChain->SetBranchStatus("eleSvCosAngle",1);
    fChain->SetBranchStatus("eleSvLxy",1);
    fChain->SetBranchStatus("eleSvLxyError",1);
    fChain->SetBranchStatus("kaonEECharge_lead",1);
    fChain->SetBranchStatus("kaonEEPt_lead",1);
    fChain->SetBranchStatus("kaonEEEta_lead",1);
    fChain->SetBranchStatus("kaonEEPhi_lead",1);
    fChain->SetBranchStatus("kaonEETrkNormChi2_lead",1);
    fChain->SetBranchStatus("kaonEECharge_sublead",1);
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
   tree->Branch("eledR",     &cutTree.eledR,     "eledR/F");
   tree->Branch("kaondR",    &cutTree.kaondR,    "kaondR/F");
   tree->Branch("jpsiPhidR", &cutTree.jpsiPhidR, "jpsiPhidR/F");
   tree->Branch("eleID_lead",&cutTree.eleID_lead,"eleID_lead/F");
   tree->Branch("eleID_sublead",&cutTree.eleID_sublead,"eleID_sublead/F");

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
      if ((jentry % 100000) == 0) cout<<"Processing entry: "<<jentry<<", accepted events: "<<nAccept<<endl;

      /***************************
       *
       * Electron Part
       *
       * ************************/

      for (int i=0; i<nEle; ++i) {

	 bool preCut = true;

	 if (eleCharge_lead->at(i) * eleCharge_sublead->at(i) < 0) continue;
	 if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) > 0) continue; // same charge
	 //if (elePt_lead->at(i) < 0.4 || elePt_sublead->at(i) < 0.4 || kaonEEPt_lead->at(i) < 0.4 || kaonEEPt_sublead->at(i) < 0.4) continue;

	 bool cutBasedSelected = true;
	 //if (elePt_lead->at(i) < 2 || elePt_sublead->at(i) < 2) cutBasedSelected = false;
	 //if (kaonEEPt_lead->at(i) < 0.8 || kaonEEPt_sublead->at(i) < 0.8) cutBasedSelected = false;
	 if (eleSvProb->at(i) < 0.01) cutBasedSelected = false;
	 //if (eleSvLxy->at(i)/eleSvLxyError->at(i) < 3) cutBasedSelected = false;
	 if (eleSvCosAngle->at(i) < 0.9) cutBasedSelected = false;
	 //if (bsEEdRele->at(i) > 1.5) cutBasedSelected = false;
	 //if (bsEEdRkaon->at(i) > 0.5) cutBasedSelected = false;
	 //if (bsEEdRJpsiPhi->at(i) > 1.5) cutBasedSelected = false;

	 double jpsiPt = TMath::Sqrt(pow(elePt_lead->at(i),2) + pow(elePt_sublead->at(i),2) + 2.0*elePt_lead->at(i)*elePt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),elePhi_sublead->at(i))))/bsEEBsMass->at(i);
	 double phiPt = TMath::Sqrt(pow(kaonEEPt_lead->at(i),2) + pow(kaonEEPt_sublead->at(i),2) + 2.0*kaonEEPt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonEEPhi_lead->at(i),kaonEEPhi_sublead->at(i))))/bsEEBsMass->at(i);
	 double bsPt = TMath::Sqrt(pow(jpsiPt*bsEEBsMass->at(i),2) + pow(phiPt*bsEEBsMass->at(i),2) 
		     + 2.0*elePt_lead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_lead->at(i))) 
		     + 2.0*elePt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_sublead->at(i))) 
		     + 2.0*elePt_sublead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_lead->at(i))) 
		     + 2.0*elePt_sublead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_sublead->at(i))) 
		     )/bsEEBsMass->at(i);

	 //if (bsPt < 3) cutBasedSelected = false;

	 if (cutBasedSelected) {

	    cutTree.m_ee = bsEEJpsiMass->at(i);
	    cutTree.m_KK = bsEEPhiMass->at(i);
	    cutTree.m_KKee = bsEEBsMass->at(i);
	    cutTree.elePtLead = elePt_lead->at(i);
	    cutTree.elePtSublead = elePt_sublead->at(i);
	    cutTree.kaonPtLead = kaonEEPt_lead->at(i);
	    cutTree.kaonPtSublead = kaonEEPt_sublead->at(i);
	    cutTree.bsPt = bsPt*bsEEBsMass->at(i);
	    cutTree.LxySig = eleSvLxy->at(i)/eleSvLxyError->at(i);
	    cutTree.SvProb = eleSvProb->at(i);
	    cutTree.SvCosine = eleSvCosAngle->at(i);
	    cutTree.eledR = bsEEdRele->at(i);
	    cutTree.kaondR = bsEEdRkaon->at(i);
	    cutTree.jpsiPhidR = bsEEdRJpsiPhi->at(i);

	    cutTree.eleID_lead = 0;
	    if (eleIDbit_lead->at(i)>>0&1) cutTree.eleID_lead = cutTree.eleID_lead + 1;
	    if (eleIDbit_lead->at(i)>>1&1) cutTree.eleID_lead = cutTree.eleID_lead + 1;
	    if (eleIDbit_lead->at(i)>>2&1) cutTree.eleID_lead = cutTree.eleID_lead + 1;
	    if (eleIDbit_lead->at(i)>>3&1) cutTree.eleID_lead = cutTree.eleID_lead + 1;

	    cutTree.eleID_sublead = 0;
	    if (eleIDbit_sublead->at(i)>>0&1) cutTree.eleID_sublead = cutTree.eleID_sublead + 1;
	    if (eleIDbit_sublead->at(i)>>1&1) cutTree.eleID_sublead = cutTree.eleID_sublead + 1;
	    if (eleIDbit_sublead->at(i)>>2&1) cutTree.eleID_sublead = cutTree.eleID_sublead + 1;
	    if (eleIDbit_sublead->at(i)>>3&1) cutTree.eleID_sublead = cutTree.eleID_sublead + 1;

	    tree->Fill();


	 }

      }

   }

   file_out.cd();
   file_out.Write();
   file_out.Close();


}
