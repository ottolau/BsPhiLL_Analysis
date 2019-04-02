#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TRandom.h"

typedef struct {
  Float_t       elePtLead;
  Float_t       eleEtaLead;
  Float_t       elePhiLead;
  Float_t       eleD0Lead;
  Float_t       eleDzLead;
  Float_t       elePtSublead;
  Float_t       eleEtaSublead;
  Float_t       elePhiSublead;
  Float_t       eleD0Sublead;
  Float_t       eleDzSublead;
  Float_t       kaonPtLead;
  Float_t       kaonEtaLead;
  Float_t       kaonPhiLead;
  Float_t       kaonD0Lead;
  Float_t       kaonDzLead;
  Float_t       kaonNormChi2Lead;
  Float_t       kaonPtSublead;
  Float_t       kaonEtaSublead;
  Float_t       kaonPhiSublead;
  Float_t       kaonD0Sublead;
  Float_t       kaonDzSublead;
  Float_t       kaonNormChi2Sublead;
  Float_t       eledR;
  Float_t       kaondR;
  Float_t       jpsiPhidR;
  Float_t       svCtxy;
  Float_t       svProb;
  Float_t       svCosine;
  Float_t       svLxy;
  Float_t       svLxySig;
  Float_t       jpsiPt;
  Float_t       phiPt;
  Float_t       bsPt;
  Float_t       jpsiMass;
  Float_t       phiMass;
  Float_t       bsMass;
  Float_t	jpsiMassFrac;
  Float_t	phiMassFrac;
  Float_t	bsMassFrac;
} bsFeatures_t;

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

   // Initialization
   TTree *tree = new TTree("background", "background tree");
   bsFeatures_t bsFeatures;
   tree->Branch("elePtLead",      &bsFeatures.elePtLead,      "elePtLead/F");
   tree->Branch("eleEtaLead",     &bsFeatures.eleEtaLead,     "eleEtaLead/F");
   tree->Branch("elePhiLead",     &bsFeatures.elePhiLead,     "elePhiLead/F");
   tree->Branch("eleD0Lead",     &bsFeatures.eleD0Lead,     "eleD0Lead/F");
   tree->Branch("eleDzLead",     &bsFeatures.eleDzLead,     "eleDzLead/F");

   tree->Branch("elePtSublead",      &bsFeatures.elePtSublead,      "elePtLead/F");
   tree->Branch("eleEtaSublead",     &bsFeatures.eleEtaSublead,     "eleEtaSublead/F");
   tree->Branch("elePhiSublead",     &bsFeatures.elePhiSublead,     "elePhiSublead/F");
   tree->Branch("eleD0Sublead",     &bsFeatures.eleD0Sublead,     "eleD0Sublead/F");
   tree->Branch("eleDzSublead",     &bsFeatures.eleDzSublead,     "eleDzSublead/F");

   tree->Branch("kaonPtLead",    &bsFeatures.kaonPtLead,    "kaonPtLead/F");
   tree->Branch("kaonEtaLead",   &bsFeatures.kaonEtaLead,   "kaonEtaLead/F");
   tree->Branch("kaonPhiLead",   &bsFeatures.kaonPhiLead,   "kaonPhiLead/F");
   tree->Branch("kaonD0Lead",   &bsFeatures.kaonD0Lead,   "kaonD0Lead/F");
   tree->Branch("kaonDzLead",   &bsFeatures.kaonDzLead,   "kaonDzLead/F");
   tree->Branch("kaonNormChi2Lead",   &bsFeatures.kaonNormChi2Lead,   "kaonNormChi2Lead/F");

   tree->Branch("kaonPtSublead",    &bsFeatures.kaonPtSublead,    "kaonPtLead/F");
   tree->Branch("kaonEtaSublead",   &bsFeatures.kaonEtaSublead,   "kaonEtaSublead/F");
   tree->Branch("kaonPhiSublead",   &bsFeatures.kaonPhiSublead,   "kaonPhiSublead/F");
   tree->Branch("kaonD0Sublead",   &bsFeatures.kaonD0Sublead,   "kaonD0Sublead/F");
   tree->Branch("kaonDzSublead",   &bsFeatures.kaonDzSublead,   "kaonDzSublead/F");
   tree->Branch("kaonNormChi2Sublead",   &bsFeatures.kaonNormChi2Sublead,   "kaonNormChi2Sublead/F");

   tree->Branch("eledR",         &bsFeatures.eledR,         "eledR/F");
   tree->Branch("kaondR",       &bsFeatures.kaondR,       "kaondR/F");
   tree->Branch("jpsiPhidR",    &bsFeatures.jpsiPhidR,    "jpsiPhidR/F");

   tree->Branch("svCtxy",       &bsFeatures.svCtxy,       "svCtxy/F");
   tree->Branch("svProb",       &bsFeatures.svProb,       "svProb/F");
   tree->Branch("svCosine",     &bsFeatures.svCosine,     "svCosine/F");
   tree->Branch("svLxy",     &bsFeatures.svLxy,     "svLxy/F");
   tree->Branch("svLxySig",     &bsFeatures.svLxySig,     "svLxySig/F");

   tree->Branch("jpsiPt",       &bsFeatures.jpsiPt,       "jpsiPt/F");
   tree->Branch("phiPt",        &bsFeatures.phiPt,        "phiPt/F");
   tree->Branch("bsPt",         &bsFeatures.bsPt,         "bsPt/F");

   tree->Branch("jpsiMass",     &bsFeatures.jpsiMass,     "jpsiMass/F");
   tree->Branch("phiMass",      &bsFeatures.phiMass,      "phiMass/F");
   tree->Branch("bsMass",       &bsFeatures.bsMass,       "bsMass/F");
   tree->Branch("jpsiMassFrac",     &bsFeatures.jpsiMassFrac,     "jpsiMassFrac/F");
   tree->Branch("phiMassFrac",      &bsFeatures.phiMassFrac,      "phiMassFrac/F");
   tree->Branch("bsMassFrac",       &bsFeatures.bsMassFrac,       "bsMassFrac/F");


   TRandom3 rand;
   rand.SetSeed(0); 
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
	 if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) < 0) continue; // same charge
	 //if (elePt_lead->at(i) < 0.4 || elePt_sublead->at(i) < 0.4 || kaonEEPt_lead->at(i) < 0.4 || kaonEEPt_sublead->at(i) < 0.4) continue;
	 if (elePt_lead->at(i) < 2.0 || elePt_sublead->at(i) < 2.0) continue;
	 if (kaonEEPt_lead->at(i) < 1.0 || kaonEEPt_sublead->at(i) < 0.7) continue;
	 if (eleSvCosAngle->at(i) < 0.9) continue;
         if (!eleConvVeto_lead->at(i) || !eleConvVeto_sublead->at(i)) continue;

	 if (preCut == true) {

    	    //if  (!((jpsisideminuslow < bsEEJpsiMass->at(i) && bsEEJpsiMass->at(i) < jpsisideplusup) && ((phisideminuslow < bsEEPhiMass->at(i) && bsEEPhiMass->at(i) < phisideminusup) || (phisidepluslow < bsEEPhiMass->at(i) && bsEEPhiMass->at(i) < phisideplusup)) && (bslow < bsEEBsMass->at(i) && bsEEBsMass->at(i) < bsup))) continue;
     	    //if  (!((philow < bsEEPhiMass->at(i) && bsEEPhiMass->at(i) < phiup) && (bslow < bsEEBsMass->at(i) && bsEEBsMass->at(i) < bsup))) continue;
  
            //if (rand.Uniform(0.0, 1.0) > 0.5) continue;

	    bsFeatures.elePtLead = elePt_lead->at(i)/bsEEBsMass->at(i);
	    bsFeatures.eleEtaLead = eleEta_lead->at(i);
	    bsFeatures.elePhiLead = elePhi_lead->at(i);
	    bsFeatures.eleD0Lead = eleD0_lead->at(i);
	    bsFeatures.eleDzLead = eleDz_lead->at(i);

	    bsFeatures.elePtSublead = elePt_sublead->at(i)/bsEEBsMass->at(i);
	    bsFeatures.eleEtaSublead = eleEta_sublead->at(i);
	    bsFeatures.elePhiSublead = elePhi_sublead->at(i);
	    bsFeatures.eleD0Sublead = eleD0_sublead->at(i);
	    bsFeatures.eleDzSublead = eleDz_sublead->at(i);

	    bsFeatures.kaonPtLead = kaonEEPt_lead->at(i)/bsEEBsMass->at(i);
	    bsFeatures.kaonEtaLead = kaonEEEta_lead->at(i);
	    bsFeatures.kaonPhiLead = kaonEEPhi_lead->at(i);
	    bsFeatures.kaonD0Lead = kaonEED0_lead->at(i);
	    bsFeatures.kaonDzLead = kaonEEDz_lead->at(i);
	    bsFeatures.kaonNormChi2Lead = kaonEETrkNormChi2_lead->at(i);

	    bsFeatures.kaonPtSublead = kaonEEPt_sublead->at(i)/bsEEBsMass->at(i);
	    bsFeatures.kaonEtaSublead = kaonEEEta_sublead->at(i);
	    bsFeatures.kaonPhiSublead = kaonEEPhi_sublead->at(i);
	    bsFeatures.kaonD0Sublead = kaonEED0_sublead->at(i);
	    bsFeatures.kaonDzSublead = kaonEEDz_sublead->at(i);
	    bsFeatures.kaonNormChi2Sublead = kaonEETrkNormChi2_sublead->at(i);

	    bsFeatures.eledR = bsEEdRele->at(i);
	    bsFeatures.kaondR = bsEEdRkaon->at(i);
	    bsFeatures.jpsiPhidR = bsEEdRJpsiPhi->at(i);

	    bsFeatures.jpsiPt = TMath::Sqrt(pow(elePt_lead->at(i),2) + pow(elePt_sublead->at(i),2) + 2.0*elePt_lead->at(i)*elePt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),elePhi_sublead->at(i))))/bsEEBsMass->at(i);
	    bsFeatures.phiPt = TMath::Sqrt(pow(kaonEEPt_lead->at(i),2) + pow(kaonEEPt_sublead->at(i),2) + 2.0*kaonEEPt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonEEPhi_lead->at(i),kaonEEPhi_sublead->at(i))))/bsEEBsMass->at(i);
	    bsFeatures.bsPt = TMath::Sqrt(pow(bsFeatures.jpsiPt*bsEEBsMass->at(i),2) + pow(bsFeatures.phiPt*bsEEBsMass->at(i),2) 
	                + 2.0*elePt_lead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_lead->at(i))) 
			+ 2.0*elePt_lead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_lead->at(i),kaonEEPhi_sublead->at(i))) 
			+ 2.0*elePt_sublead->at(i)*kaonEEPt_lead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_lead->at(i))) 
			+ 2.0*elePt_sublead->at(i)*kaonEEPt_sublead->at(i)*TMath::Cos(deltaPhi(elePhi_sublead->at(i),kaonEEPhi_sublead->at(i))) 
			)/bsEEBsMass->at(i);
	    bsFeatures.svCtxy = eleSvCtxy->at(i);
	    bsFeatures.svProb = eleSvProb->at(i);
	    bsFeatures.svCosine = eleSvCosAngle->at(i);
	    bsFeatures.svLxy = eleSvLxy->at(i);
	    bsFeatures.svLxySig = eleSvLxy->at(i)/eleSvLxyError->at(i);

	    bsFeatures.jpsiMass = bsEEJpsiMass->at(i);
	    bsFeatures.phiMass = bsEEPhiMass->at(i);
	    bsFeatures.bsMass = bsEEBsMass->at(i);
	    bsFeatures.jpsiMassFrac = bsEEJpsiMass->at(i)/jpsiM;
	    bsFeatures.phiMassFrac = bsEEPhiMass->at(i)/phiM;
	    bsFeatures.bsMassFrac = bsEEBsMass->at(i)/bsM;

	    tree->Fill();
	    nAccept++;
	    break;
	 }
      }

   }

   file_out.cd();
   file_out.Write();
   file_out.Close();


}
