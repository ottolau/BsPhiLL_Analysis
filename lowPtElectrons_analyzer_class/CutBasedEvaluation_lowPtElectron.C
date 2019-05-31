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

   // Histograms
   // Electron
 
   TH1D *h_q2_InvM_pass = new TH1D("h_q2_InvM_pass", "", 200, 0.0, 25.0);
   TH1D *h_jpsi_InvM_pass = new TH1D("h_jpsi_InvM_pass", "", 100, 2.4, 3.8);
   TH1D *h_phiee_InvM_pass = new TH1D("h_phiee_InvM_pass", "", 200, 0.98, 1.1);
   TH1D *h_bsee_InvM_pass = new TH1D("h_bsee_InvM_pass", "", 200, 4.5, 6.0);
  
   TH1D *h_electron_pt_lead_pass = new TH1D("h_electron_pt_lead_pass", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_lead_pass = new TH1D("h_electron_eta_lead_pass", "", 50, -3, 3);
   TH1D *h_electron_phi_lead_pass = new TH1D("h_electron_phi_lead_pass", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_lead_pass = new TH1D("h_electron_D0_lead_pass", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_lead_pass = new TH1D("h_electron_Dz_lead_pass", "", 100, -1.0, 1.0);
   TH1D *h_electron_MVABWP_lead_pass = new TH1D("h_electron_MVABWP_lead_pass", "", 100, -10.0, 10.0);
   TH1D *h_electron_MVAUnBWP_lead_pass = new TH1D("h_electron_MVAUnBWP_lead_pass", "", 100, -10.0, 10.0);

   TH1D *h_electron_pt_sublead_pass = new TH1D("h_electron_pt_sublead_pass", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_sublead_pass = new TH1D("h_electron_eta_sublead_pass", "", 50, -3, 3);
   TH1D *h_electron_phi_sublead_pass = new TH1D("h_electron_phi_sublead_pass", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_sublead_pass = new TH1D("h_electron_D0_sublead_pass", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_sublead_pass = new TH1D("h_electron_Dz_sublead_pass", "", 100, -1.0, 1.0);
   TH1D *h_electron_MVABWP_sublead_pass = new TH1D("h_electron_MVABWP_sublead_pass", "", 100, -10.0, 10.0);
   TH1D *h_electron_MVAUnBWP_sublead_pass = new TH1D("h_electron_MVAUnBWP_sublead_pass", "", 100, -10.0, 10.0);

   TH1D *h_kaon_ee_pt_lead_pass = new TH1D("h_kaon_ee_pt_lead_pass", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_lead_pass = new TH1D("h_kaon_ee_eta_lead_pass", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_lead_pass = new TH1D("h_kaon_ee_phi_lead_pass", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_lead_pass = new TH1D("h_kaon_ee_Dz_lead_pass", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_lead_pass = new TH1D("h_kaon_ee_D0_lead_pass", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_lead_pass = new TH1D("h_kaon_ee_D0Sig_lead_pass", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_lead_pass = new TH1D("h_kaon_ee_normChi2_lead_pass", "", 100, 0.0, 5.0);

   TH1D *h_kaon_ee_pt_sublead_pass = new TH1D("h_kaon_ee_pt_sublead_pass", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_sublead_pass = new TH1D("h_kaon_ee_eta_sublead_pass", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_sublead_pass = new TH1D("h_kaon_ee_phi_sublead_pass", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_sublead_pass = new TH1D("h_kaon_ee_Dz_sublead_pass", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_sublead_pass = new TH1D("h_kaon_ee_D0_sublead_pass", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_sublead_pass = new TH1D("h_kaon_ee_D0Sig_sublead_pass", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_sublead_pass = new TH1D("h_kaon_ee_normChi2_sublead_pass", "", 100, 0.0, 5.0);

   TH1D *h_bs_pt_ee_pass = new TH1D("h_bs_pt_ee_pass", "", 100, 0.0, 50.0);

   TH1D *h_svProb_ee_pass = new TH1D("h_svProb_ee_pass", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_ee_pass = new TH1D("h_svCosAngle_ee_pass", "", 200, -1.0, 1.0);
   TH1D *h_svCtxy_ee_pass = new TH1D("h_svCtxy_ee_pass", "", 100, -1.0, 1.0);
   TH1D *h_svLxy_ee_pass = new TH1D("h_svLxy_ee_pass", "", 100, 0.0, 1);
   TH1D *h_svLxySig_ee_pass = new TH1D("h_svLxySig_ee_pass", "", 200, 0.0, 60);

   TH1D *h_electron_dR_pass = new TH1D("h_electron_dR_pass", "", 120, 0.0, 4.0);
   TH1D *h_kaon_ee_dR_pass = new TH1D("h_kaon_ee_dR_pass", "", 120, 0.0, 4.0);
   TH1D *h_jpsiPhiOpen_ee_pass = new TH1D("h_jpsiPhiOpen_ee_pass", "", 120, 0.0, 4.0);


   TH1D *h_q2_InvM_veto = new TH1D("h_q2_InvM_veto", "", 200, 0.0, 25.0);
   TH1D *h_jpsi_InvM_veto = new TH1D("h_jpsi_InvM_veto", "", 100, 2.4, 3.8);
   TH1D *h_phiee_InvM_veto = new TH1D("h_phiee_InvM_veto", "", 200, 0.98, 1.1);
   TH1D *h_bsee_InvM_veto = new TH1D("h_bsee_InvM_veto", "", 200, 4.5, 6.0);
  
   TH1D *h_electron_pt_lead_veto = new TH1D("h_electron_pt_lead_veto", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_lead_veto = new TH1D("h_electron_eta_lead_veto", "", 50, -3, 3);
   TH1D *h_electron_phi_lead_veto = new TH1D("h_electron_phi_lead_veto", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_lead_veto = new TH1D("h_electron_D0_lead_veto", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_lead_veto = new TH1D("h_electron_Dz_lead_veto", "", 100, -1.0, 1.0);
   TH1D *h_electron_MVABWP_lead_veto = new TH1D("h_electron_MVABWP_lead_veto", "", 100, -10.0, 10.0);
   TH1D *h_electron_MVAUnBWP_lead_veto = new TH1D("h_electron_MVAUnBWP_lead_veto", "", 100, -10.0, 10.0);

   TH1D *h_electron_pt_sublead_veto = new TH1D("h_electron_pt_sublead_veto", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_sublead_veto = new TH1D("h_electron_eta_sublead_veto", "", 50, -3, 3);
   TH1D *h_electron_phi_sublead_veto = new TH1D("h_electron_phi_sublead_veto", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_sublead_veto = new TH1D("h_electron_D0_sublead_veto", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_sublead_veto = new TH1D("h_electron_Dz_sublead_veto", "", 100, -1.0, 1.0);
   TH1D *h_electron_MVABWP_sublead_veto = new TH1D("h_electron_MVABWP_sublead_veto", "", 100, -10.0, 10.0);
   TH1D *h_electron_MVAUnBWP_sublead_veto = new TH1D("h_electron_MVAUnBWP_sublead_veto", "", 100, -10.0, 10.0);

   TH1D *h_kaon_ee_pt_lead_veto = new TH1D("h_kaon_ee_pt_lead_veto", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_lead_veto = new TH1D("h_kaon_ee_eta_lead_veto", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_lead_veto = new TH1D("h_kaon_ee_phi_lead_veto", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_lead_veto = new TH1D("h_kaon_ee_Dz_lead_veto", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_lead_veto = new TH1D("h_kaon_ee_D0_lead_veto", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_lead_veto = new TH1D("h_kaon_ee_D0Sig_lead_veto", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_lead_veto = new TH1D("h_kaon_ee_normChi2_lead_veto", "", 100, 0.0, 5.0);

   TH1D *h_kaon_ee_pt_sublead_veto = new TH1D("h_kaon_ee_pt_sublead_veto", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_sublead_veto = new TH1D("h_kaon_ee_eta_sublead_veto", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_sublead_veto = new TH1D("h_kaon_ee_phi_sublead_veto", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_sublead_veto = new TH1D("h_kaon_ee_Dz_sublead_veto", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_sublead_veto = new TH1D("h_kaon_ee_D0_sublead_veto", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_sublead_veto = new TH1D("h_kaon_ee_D0Sig_sublead_veto", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_sublead_veto = new TH1D("h_kaon_ee_normChi2_sublead_veto", "", 100, 0.0, 5.0);

   TH1D *h_bs_pt_ee_veto = new TH1D("h_bs_pt_ee_veto", "", 100, 0.0, 50.0);

   TH1D *h_svProb_ee_veto = new TH1D("h_svProb_ee_veto", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_ee_veto = new TH1D("h_svCosAngle_ee_veto", "", 200, -1.0, 1.0);
   TH1D *h_svCtxy_ee_veto = new TH1D("h_svCtxy_ee_veto", "", 100, -1.0, 1.0);
   TH1D *h_svLxy_ee_veto = new TH1D("h_svLxy_ee_veto", "", 100, 0.0, 1);
   TH1D *h_svLxySig_ee_veto = new TH1D("h_svLxySig_ee_veto", "", 200, 0.0, 60);

   TH1D *h_electron_dR_veto = new TH1D("h_electron_dR_veto", "", 120, 0.0, 4.0);
   TH1D *h_kaon_ee_dR_veto = new TH1D("h_kaon_ee_dR_veto", "", 120, 0.0, 4.0);
   TH1D *h_jpsiPhiOpen_ee_veto = new TH1D("h_jpsiPhiOpen_ee_veto", "", 120, 0.0, 4.0);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if ((jentry%10000) == 0) std::cout<<"Procssing entry: "<<jentry<<std::endl;

      vector< pair<float,float> > jpsiPtList_ee_pass;
      vector< pair<float,float> > phiPtList_ee_pass;
      vector< pair<float,float> > jpsiPtList_ee_veto;
      vector< pair<float,float> > phiPtList_ee_veto;

      /***************************
       *
       * Electron Part
       *
       * ************************/

      for (int i=0; i<nLowPt; ++i) {

	 bool preCut = true;

	 if (lowPtCharge_lead->at(i) * lowPtCharge_sublead->at(i) > 0) continue;
	 if (kaonLowPtCharge_lead->at(i) * kaonLowPtCharge_sublead->at(i) > 0) continue; // same charge
	 //if (lowPtPt_lead->at(i) < 0.4 || lowPtPt_sublead->at(i) < 0.4 || kaonLowPtPt_lead->at(i) < 0.4 || kaonLowPtPt_sublead->at(i) < 0.4) continue;

	 pair<float,float> jpsiPt_ee = make_pair(lowPtPt_lead->at(i), lowPtPt_sublead->at(i));
	 pair<float,float> phiPt_ee = make_pair(kaonLowPtPt_lead->at(i), kaonLowPtPt_sublead->at(i));

	 bool cutBasedSelected = true;
	 //if (lowPtPt_lead->at(i) < 2 || lowPtPt_sublead->at(i) < 2) cutBasedSelected = false;
	 if (kaonLowPtPt_lead->at(i) < 0.8 || kaonLowPtPt_sublead->at(i) < 0.4) cutBasedSelected = false;
	 if (lowPtSvProb->at(i) < 0.01) cutBasedSelected = false;
	 if (lowPtSvLxy->at(i)/lowPtSvLxyError->at(i) < 1) cutBasedSelected = false;
	 if (lowPtSvCosAngle->at(i) < 0.9) cutBasedSelected = false;
	 //if (bsLowPtdRele->at(i) > 1.5) cutBasedSelected = false;
	 //if (bsLowPtdRkaon->at(i) > 0.5) cutBasedSelected = false;
	 //if (bsLowPtdRJpsiPhi->at(i) > 1.5) cutBasedSelected = false;

	 double jpsiPt = TMath::Sqrt(pow(lowPtPt_lead->at(i),2) + pow(lowPtPt_sublead->at(i),2) + 2.0*lowPtPt_lead->at(i)*lowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),lowPtPhi_sublead->at(i))))/bsLowPtBsMass->at(i);
	 double phiPt = TMath::Sqrt(pow(kaonLowPtPt_lead->at(i),2) + pow(kaonLowPtPt_sublead->at(i),2) + 2.0*kaonLowPtPt_lead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(kaonLowPtPhi_lead->at(i),kaonLowPtPhi_sublead->at(i))))/bsLowPtBsMass->at(i);
	 double bsPt = TMath::Sqrt(pow(jpsiPt*bsLowPtBsMass->at(i),2) + pow(phiPt*bsLowPtBsMass->at(i),2) 
		     + 2.0*lowPtPt_lead->at(i)*kaonLowPtPt_lead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),kaonLowPtPhi_lead->at(i))) 
		     + 2.0*lowPtPt_lead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_lead->at(i),kaonLowPtPhi_sublead->at(i))) 
		     + 2.0*lowPtPt_sublead->at(i)*kaonLowPtPt_lead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_sublead->at(i),kaonLowPtPhi_lead->at(i))) 
		     + 2.0*lowPtPt_sublead->at(i)*kaonLowPtPt_sublead->at(i)*TMath::Cos(deltaPhi(lowPtPhi_sublead->at(i),kaonLowPtPhi_sublead->at(i))) 
		     )/bsLowPtBsMass->at(i);

	 if (bsPt < 6) cutBasedSelected = false;

	 if (cutBasedSelected) {

	    if (find(jpsiPtList_ee_pass.begin(), jpsiPtList_ee_pass.end(), jpsiPt_ee) == jpsiPtList_ee_pass.end()) {
	       //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	       //}
	       h_electron_pt_lead_pass->Fill(lowPtPt_lead->at(i));
	       h_electron_eta_lead_pass->Fill(lowPtEta_lead->at(i));
	       h_electron_phi_lead_pass->Fill(lowPtPhi_lead->at(i));
	       h_electron_D0_lead_pass->Fill(lowPtD0_lead->at(i));
	       h_electron_Dz_lead_pass->Fill(lowPtDz_lead->at(i));
	       h_electron_MVABWP_lead_pass->Fill(lowPtMVABWP_lead->at(i));
	       h_electron_MVAUnBWP_lead_pass->Fill(lowPtMVAUnBWP_lead->at(i));

	       h_electron_pt_sublead_pass->Fill(lowPtPt_sublead->at(i));
	       h_electron_eta_sublead_pass->Fill(lowPtEta_sublead->at(i));
	       h_electron_phi_sublead_pass->Fill(lowPtPhi_sublead->at(i));
	       h_electron_D0_sublead_pass->Fill(lowPtD0_sublead->at(i));
	       h_electron_Dz_sublead_pass->Fill(lowPtDz_sublead->at(i));
	       h_electron_MVABWP_sublead_pass->Fill(lowPtMVABWP_sublead->at(i));
	       h_electron_MVAUnBWP_sublead_pass->Fill(lowPtMVAUnBWP_sublead->at(i));

	       h_electron_dR_pass->Fill(bsLowPtdRele->at(i));
	       h_q2_InvM_pass->Fill(bsLowPtJpsiMass->at(i)*bsLowPtJpsiMass->at(i));
	       h_jpsi_InvM_pass->Fill(bsLowPtJpsiMass->at(i));

	       jpsiPtList_ee_pass.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee_pass.begin(), phiPtList_ee_pass.end(), phiPt_ee) == phiPtList_ee_pass.end()) {
	       //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	       //}
	       h_kaon_ee_pt_lead_pass->Fill(kaonLowPtPt_lead->at(i));
	       h_kaon_ee_eta_lead_pass->Fill(kaonLowPtEta_lead->at(i));
	       h_kaon_ee_phi_lead_pass->Fill(kaonLowPtPhi_lead->at(i));
	       h_kaon_ee_Dz_lead_pass->Fill(kaonLowPtDz_lead->at(i));
	       h_kaon_ee_D0_lead_pass->Fill(kaonLowPtD0_lead->at(i));
	       h_kaon_ee_D0Sig_lead_pass->Fill(kaonLowPtD0_lead->at(i)/kaonLowPtD0Error_lead->at(i));
	       h_kaon_ee_normChi2_lead_pass->Fill(kaonLowPtTrkNormChi2_lead->at(i));

	       h_kaon_ee_pt_sublead_pass->Fill(kaonLowPtPt_sublead->at(i));
	       h_kaon_ee_eta_sublead_pass->Fill(kaonLowPtEta_sublead->at(i));
	       h_kaon_ee_phi_sublead_pass->Fill(kaonLowPtPhi_sublead->at(i));
	       h_kaon_ee_Dz_sublead_pass->Fill(kaonLowPtDz_sublead->at(i));
	       h_kaon_ee_D0_sublead_pass->Fill(kaonLowPtD0_sublead->at(i));
	       h_kaon_ee_D0Sig_sublead_pass->Fill(kaonLowPtD0_sublead->at(i)/kaonLowPtD0Error_sublead->at(i));
	       h_kaon_ee_normChi2_sublead_pass->Fill(kaonLowPtTrkNormChi2_sublead->at(i));

	       h_kaon_ee_dR_pass->Fill(bsLowPtdRkaon->at(i));
	       h_phiee_InvM_pass->Fill(bsLowPtPhiMass->at(i));

	       phiPtList_ee_pass.push_back(phiPt_ee);
	    }
	    //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	    //}
	    
	    h_bs_pt_ee_pass->Fill(bsPt);
	    h_svProb_ee_pass->Fill(lowPtSvProb->at(i));
	    h_svCosAngle_ee_pass->Fill(lowPtSvCosAngle->at(i));
	    h_svCtxy_ee_pass->Fill(lowPtSvCtxy->at(i));
	    h_svLxy_ee_pass->Fill(lowPtSvLxy->at(i));
	    h_svLxySig_ee_pass->Fill(lowPtSvLxy->at(i)/lowPtSvLxyError->at(i));

	    h_jpsiPhiOpen_ee_pass->Fill(bsLowPtdRJpsiPhi->at(i));

	    if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup) {
	       h_bsee_InvM_pass->Fill(bsLowPtBsMass->at(i));
	    }

	 } else {

	    if (find(jpsiPtList_ee_veto.begin(), jpsiPtList_ee_veto.end(), jpsiPt_ee) == jpsiPtList_ee_veto.end()) {
	       //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	       //}
	       h_electron_pt_lead_veto->Fill(lowPtPt_lead->at(i));
	       h_electron_eta_lead_veto->Fill(lowPtEta_lead->at(i));
	       h_electron_phi_lead_veto->Fill(lowPtPhi_lead->at(i));
	       h_electron_D0_lead_veto->Fill(lowPtD0_lead->at(i));
	       h_electron_Dz_lead_veto->Fill(lowPtDz_lead->at(i));
	       h_electron_MVABWP_lead_veto->Fill(lowPtMVABWP_lead->at(i));
	       h_electron_MVAUnBWP_lead_veto->Fill(lowPtMVAUnBWP_lead->at(i));

	       h_electron_pt_sublead_veto->Fill(lowPtPt_sublead->at(i));
	       h_electron_eta_sublead_veto->Fill(lowPtEta_sublead->at(i));
	       h_electron_phi_sublead_veto->Fill(lowPtPhi_sublead->at(i));
	       h_electron_D0_sublead_veto->Fill(lowPtD0_sublead->at(i));
	       h_electron_Dz_sublead_veto->Fill(lowPtDz_sublead->at(i));
	       h_electron_MVABWP_sublead_veto->Fill(lowPtMVABWP_sublead->at(i));
	       h_electron_MVAUnBWP_sublead_veto->Fill(lowPtMVAUnBWP_sublead->at(i));

	       h_electron_dR_veto->Fill(bsLowPtdRele->at(i));
	       h_q2_InvM_veto->Fill(bsLowPtJpsiMass->at(i)*bsLowPtJpsiMass->at(i));
	       h_jpsi_InvM_veto->Fill(bsLowPtJpsiMass->at(i));

	       jpsiPtList_ee_veto.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee_veto.begin(), phiPtList_ee_veto.end(), phiPt_ee) == phiPtList_ee_veto.end()) {
	       //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	       //}
	       h_kaon_ee_pt_lead_veto->Fill(kaonLowPtPt_lead->at(i));
	       h_kaon_ee_eta_lead_veto->Fill(kaonLowPtEta_lead->at(i));
	       h_kaon_ee_phi_lead_veto->Fill(kaonLowPtPhi_lead->at(i));
	       h_kaon_ee_Dz_lead_veto->Fill(kaonLowPtDz_lead->at(i));
	       h_kaon_ee_D0_lead_veto->Fill(kaonLowPtD0_lead->at(i));
	       h_kaon_ee_D0Sig_lead_veto->Fill(kaonLowPtD0_lead->at(i)/kaonLowPtD0Error_lead->at(i));
	       h_kaon_ee_normChi2_lead_veto->Fill(kaonLowPtTrkNormChi2_lead->at(i));

	       h_kaon_ee_pt_sublead_veto->Fill(kaonLowPtPt_sublead->at(i));
	       h_kaon_ee_eta_sublead_veto->Fill(kaonLowPtEta_sublead->at(i));
	       h_kaon_ee_phi_sublead_veto->Fill(kaonLowPtPhi_sublead->at(i));
	       h_kaon_ee_Dz_sublead_veto->Fill(kaonLowPtDz_sublead->at(i));
	       h_kaon_ee_D0_sublead_veto->Fill(kaonLowPtD0_sublead->at(i));
	       h_kaon_ee_D0Sig_sublead_veto->Fill(kaonLowPtD0_sublead->at(i)/kaonLowPtD0Error_sublead->at(i));
	       h_kaon_ee_normChi2_sublead_veto->Fill(kaonLowPtTrkNormChi2_sublead->at(i));

	       h_kaon_ee_dR_veto->Fill(bsLowPtdRkaon->at(i));
	       h_phiee_InvM_veto->Fill(bsLowPtPhiMass->at(i));

	       phiPtList_ee_veto.push_back(phiPt_ee);
	    }
	    //if (bsLowPtJpsiMass->at(i) > jpsilow && bsLowPtJpsiMass->at(i) < jpsiup && bsLowPtPhiMass->at(i) > philow && bsLowPtPhiMass->at(i) < phiup && bsLowPtBsMass->at(i) > bslow && bsLowPtBsMass->at(i) < bsup) {
	    //}

	    h_bs_pt_ee_veto->Fill(bsPt);
	    h_svProb_ee_veto->Fill(lowPtSvProb->at(i));
	    h_svCosAngle_ee_veto->Fill(lowPtSvCosAngle->at(i));
	    h_svCtxy_ee_veto->Fill(lowPtSvCtxy->at(i));
	    h_svLxy_ee_veto->Fill(lowPtSvLxy->at(i));
	    h_svLxySig_ee_veto->Fill(lowPtSvLxy->at(i)/lowPtSvLxyError->at(i));

	    h_jpsiPhiOpen_ee_veto->Fill(bsLowPtdRJpsiPhi->at(i));
	    h_bsee_InvM_veto->Fill(bsLowPtBsMass->at(i));
	    

	 }

      }


   }

   file_out.cd();
   file_out.Write();
   file_out.Close();

}
