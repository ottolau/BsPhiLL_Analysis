#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TRandom.h"

typedef struct {
  Float_t       muPt1st;
  Float_t       muEta1st;
  Float_t       muPhi1st;
  Float_t       muPt2nd;
  Float_t       muEta2nd;
  Float_t       muPhi2nd;
  Float_t       kaonPt1st;
  Float_t       kaonEta1st;
  Float_t       kaonPhi1st;
  Float_t       kaonPt2nd;
  Float_t       kaonEta2nd;
  Float_t       kaonPhi2nd;
  Float_t       mudR;
  Float_t       kaondR;
  Float_t       jpsiPhidR;
  Float_t       svCtxy;
  Float_t       svChi2;
  Float_t       svCosine;
  Float_t       jpsiPt;
  Float_t       phiPt;
  Float_t       bsPt;
  Float_t       jpsiMass;
  Float_t       phiMass;
  Float_t       bsMass;
  Float_t       muMva1st;
  Float_t       muMva2nd;
  Float_t       muSoftMva1st;
  Float_t       muSoftMva2nd;
  Float_t       muIDCutBased1st;
  Float_t       muIDCutBased2nd;
} bsFeatures_t;

void BsPhiLLTupleTree::Loop(TString outputfile, Int_t maxevents, Float_t mvaCut)
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
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile file_out(outputfile, "recreate");

   // Initialization
   TTree *tree = new TTree("background", "background tree");
   bsFeatures_t bsFeatures;
   tree->Branch("muPt1st",      &bsFeatures.muPt1st,      "muPt1st/F");
   tree->Branch("muEta1st",     &bsFeatures.muEta1st,     "muEta1st/F");
   tree->Branch("muPhi1st",     &bsFeatures.muPhi1st,     "muPhi1st/F");
   tree->Branch("muPt2nd",      &bsFeatures.muPt2nd,      "muPt1st/F");
   tree->Branch("muEta2nd",     &bsFeatures.muEta2nd,     "muEta2nd/F");
   tree->Branch("muPhi2nd",     &bsFeatures.muPhi2nd,     "muPhi2nd/F");

   tree->Branch("kaonPt1st",    &bsFeatures.kaonPt1st,    "kaonPt1st/F");
   tree->Branch("kaonEta1st",   &bsFeatures.kaonEta1st,   "kaonEta1st/F");
   tree->Branch("kaonPhi1st",   &bsFeatures.kaonPhi1st,   "kaonPhi1st/F");
   tree->Branch("kaonPt2nd",    &bsFeatures.kaonPt2nd,    "kaonPt1st/F");
   tree->Branch("kaonEta2nd",   &bsFeatures.kaonEta2nd,   "kaonEta2nd/F");
   tree->Branch("kaonPhi2nd",   &bsFeatures.kaonPhi2nd,   "kaonPhi2nd/F");

   tree->Branch("mudR",         &bsFeatures.mudR,         "mudR/F");
   tree->Branch("kaondR",       &bsFeatures.kaondR,       "kaondR/F");
   tree->Branch("jpsiPhidR",    &bsFeatures.jpsiPhidR,    "jpsiPhidR/F");

   tree->Branch("svCtxy",       &bsFeatures.svCtxy,       "svCtxy/F");
   tree->Branch("svChi2",       &bsFeatures.svChi2,       "svChi2/F");
   tree->Branch("svCosine",     &bsFeatures.svCosine,     "svCosine/F");

   tree->Branch("jpsiPt",       &bsFeatures.jpsiPt,       "jpsiPt/F");
   tree->Branch("phiPt",        &bsFeatures.phiPt,        "phiPt/F");
   tree->Branch("bsPt",         &bsFeatures.bsPt,         "bsPt/F");

   tree->Branch("jpsiMass",     &bsFeatures.jpsiMass,     "jpsiMass/F");
   tree->Branch("phiMass",      &bsFeatures.phiMass,      "phiMass/F");
   tree->Branch("bsMass",       &bsFeatures.bsMass,       "bsMass/F");

   tree->Branch("muMva1st",     &bsFeatures.muMva1st,     "muMva1st/F");
   tree->Branch("muMva2nd",     &bsFeatures.muMva2nd,     "muMva2nd/F");
   tree->Branch("muSoftMva1st", &bsFeatures.muSoftMva1st, "muSoftMva1st/F");
   tree->Branch("muSoftMva2nd", &bsFeatures.muSoftMva2nd, "muSoftMva2nd/F");
   tree->Branch("muIDCutBased1st",     &bsFeatures.muIDCutBased1st,     "muIDCutBased1st/F");
   tree->Branch("muIDCutBased2nd",     &bsFeatures.muIDCutBased2nd,     "muIDCutBased2nd/F");


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
       * Muon Part
       *
       * ************************/

      for (int i=0; i<nMu; ++i) {


	 // Non-trigger muon part
	 bool preCut_trgmu = true;

	 if (preCut_trgmu == true && (muFiredTrgs_lead->at(i) == true || muFiredTrgs_sublead->at(i) == true)) {


    	    if  (!((jpsisideminuslow < bsMMJpsiMass->at(i) && bsMMJpsiMass->at(i) < jpsisideplusup) && ((phisideminuslow < bsMMPhiMass->at(i) && bsMMPhiMass->at(i) < phisideminusup) || (phisidepluslow < bsMMPhiMass->at(i) && bsMMPhiMass->at(i) < phisideplusup)) && (bslow < bsMMBsMass->at(i) && bsMMBsMass->at(i) < bsup))) continue;
   
            //if (rand.Uniform(0.0, 1.0) > 0.5) continue;


	    bsFeatures.muPt1st = muPt_lead->at(i)/bsMMBsMass->at(i);
	    bsFeatures.muEta1st = muEta_lead->at(i);
	    bsFeatures.muPhi1st = muPhi_lead->at(i);
	    bsFeatures.muPt2nd = muPt_sublead->at(i)/bsMMBsMass->at(i);
	    bsFeatures.muEta2nd = muEta_sublead->at(i);
	    bsFeatures.muPhi2nd = muPhi_sublead->at(i);

	    bsFeatures.kaonPt1st = kaonMMPt_lead->at(i)/bsMMBsMass->at(i);
	    bsFeatures.kaonEta1st = kaonMMEta_lead->at(i);
	    bsFeatures.kaonPhi1st = kaonMMPhi_lead->at(i);
	    bsFeatures.kaonPt2nd = kaonMMPt_sublead->at(i)/bsMMBsMass->at(i);
	    bsFeatures.kaonEta2nd = kaonMMEta_sublead->at(i);
	    bsFeatures.kaonPhi2nd = kaonMMPhi_sublead->at(i);
	    bsFeatures.mudR = bsMMdRmu->at(i);
	    bsFeatures.kaondR = bsMMdRkaon->at(i);
	    bsFeatures.jpsiPhidR = bsMMdRJpsiPhi->at(i);

	    bsFeatures.jpsiPt = TMath::Sqrt(pow(muPt_lead->at(i),2) + pow(muPt_sublead->at(i),2) + 2.0*muPt_lead->at(i)*muPt_sublead->at(i)*TMath::Cos(deltaPhi(muPhi_lead->at(i),muPhi_sublead->at(i))))/bsMMBsMass->at(i);
	    bsFeatures.phiPt = TMath::Sqrt(pow(kaonMMPt_lead->at(i),2) + pow(kaonMMPt_sublead->at(i),2) + 2.0*kaonMMPt_lead->at(i)*kaonMMPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonMMPhi_lead->at(i),kaonMMPhi_sublead->at(i))))/bsMMBsMass->at(i);
	    bsFeatures.bsPt = TMath::Sqrt(pow(bsFeatures.jpsiPt*bsMMBsMass->at(i),2) + pow(bsFeatures.phiPt*bsMMBsMass->at(i),2) 
	                + 2.0*muPt_lead->at(i)*kaonMMPt_lead->at(i)*TMath::Cos(deltaPhi(muPhi_lead->at(i),kaonMMPhi_lead->at(i))) 
			+ 2.0*muPt_lead->at(i)*kaonMMPt_sublead->at(i)*TMath::Cos(deltaPhi(muPhi_lead->at(i),kaonMMPhi_sublead->at(i))) 
			+ 2.0*muPt_sublead->at(i)*kaonMMPt_lead->at(i)*TMath::Cos(deltaPhi(muPhi_sublead->at(i),kaonMMPhi_lead->at(i))) 
			+ 2.0*muPt_sublead->at(i)*kaonMMPt_sublead->at(i)*TMath::Cos(deltaPhi(muPhi_sublead->at(i),kaonMMPhi_sublead->at(i))) 
			)/bsMMBsMass->at(i);
	    bsFeatures.svCtxy = muSvCtxy->at(i);
	    bsFeatures.svChi2 = muSvChi2->at(i)/muSvNDOF->at(i);
	    bsFeatures.svCosine = muSvCosAngle->at(i);

	    bsFeatures.jpsiMass = bsMMJpsiMass->at(i);
	    bsFeatures.phiMass = bsMMPhiMass->at(i);
	    bsFeatures.bsMass = bsMMBsMass->at(i);

	    bsFeatures.muMva1st = muIDPatMVA_lead->at(i);
	    bsFeatures.muMva2nd = muIDPatMVA_sublead->at(i);
	    bsFeatures.muSoftMva1st = muIDPatSoftMVA_lead->at(i);
	    bsFeatures.muSoftMva2nd = muIDPatSoftMVA_sublead->at(i);

	    bsFeatures.muIDCutBased1st = 0;
	    if (muIDbit_lead->at(i)>>0&1) bsFeatures.muIDCutBased1st = bsFeatures.muIDCutBased1st + 1;
	    if (muIDbit_lead->at(i)>>1&1) bsFeatures.muIDCutBased1st = bsFeatures.muIDCutBased1st + 1;
	    if (muIDbit_lead->at(i)>>2&1) bsFeatures.muIDCutBased1st = bsFeatures.muIDCutBased1st + 1;
	    bsFeatures.muIDCutBased2nd = 0;
	    if (muIDbit_sublead->at(i)>>0&1) bsFeatures.muIDCutBased2nd = bsFeatures.muIDCutBased2nd + 1;
	    if (muIDbit_sublead->at(i)>>1&1) bsFeatures.muIDCutBased2nd = bsFeatures.muIDCutBased2nd + 1;
	    if (muIDbit_sublead->at(i)>>2&1) bsFeatures.muIDCutBased2nd = bsFeatures.muIDCutBased2nd + 1;


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
