#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

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

   // Histograms
   // Electron
 
   TH1D *h_q2_InvM = new TH1D("h_q2_InvM", "", 200, 0.0, 25.0);
   TH1D *h_phiee_InvM = new TH1D("h_phiee_InvM", "", 200, 0.98, 1.1);
   TH1D *h_bsee_InvM = new TH1D("h_bsee_InvM", "", 200, 4.5, 6.0);
  
   TH1D *h_electron_pt_lead = new TH1D("h_electron_pt_lead", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_lead = new TH1D("h_electron_eta_lead", "", 50, -3, 3);
   TH1D *h_electron_phi_lead = new TH1D("h_electron_phi_lead", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_lead = new TH1D("h_electron_D0_lead", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_lead = new TH1D("h_electron_Dz_lead", "", 100, -1.0, 1.0);

   TH1D *h_electron_pt_sublead = new TH1D("h_electron_pt_sublead", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_sublead = new TH1D("h_electron_eta_sublead", "", 50, -3, 3);
   TH1D *h_electron_phi_sublead = new TH1D("h_electron_phi_sublead", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_sublead = new TH1D("h_electron_D0_sublead", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_sublead = new TH1D("h_electron_Dz_sublead", "", 100, -1.0, 1.0);

   TH1D *h_kaon_ee_pt_lead = new TH1D("h_kaon_ee_pt_lead", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_lead = new TH1D("h_kaon_ee_eta_lead", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_lead = new TH1D("h_kaon_ee_phi_lead", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_lead = new TH1D("h_kaon_ee_Dz_lead", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_lead = new TH1D("h_kaon_ee_D0_lead", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_lead = new TH1D("h_kaon_ee_D0Sig_lead", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_lead = new TH1D("h_kaon_ee_normChi2_lead", "", 100, 0.0, 5.0);

   TH1D *h_kaon_ee_pt_sublead = new TH1D("h_kaon_ee_pt_sublead", "", 80, 0.0, 5.0);
   TH1D *h_kaon_ee_eta_sublead = new TH1D("h_kaon_ee_eta_sublead", "", 50, -3, 3);
   TH1D *h_kaon_ee_phi_sublead = new TH1D("h_kaon_ee_phi_sublead", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_ee_Dz_sublead = new TH1D("h_kaon_ee_Dz_sublead", "", 100, -1.0, 1.0);
   TH1D *h_kaon_ee_D0_sublead = new TH1D("h_kaon_ee_D0_sublead", "", 100, -0.5, 0.5);
   TH1D *h_kaon_ee_D0Sig_sublead = new TH1D("h_kaon_ee_D0Sig_sublead", "", 100, -20.0, 20.0);
   TH1D *h_kaon_ee_normChi2_sublead = new TH1D("h_kaon_ee_normChi2_sublead", "", 100, 0.0, 5.0);

   TH1D *h_svProb_ee = new TH1D("h_svProb_ee", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_ee = new TH1D("h_svCosAngle_ee", "", 20, -1.0, 1.0);
   TH1D *h_svCtxy_ee = new TH1D("h_svCtxy_ee", "", 100, -1.0, 1.0);
   TH1D *h_svLxy_ee = new TH1D("h_svLxy_ee", "", 100, 0.0, 1);
   TH1D *h_svLxySig_ee = new TH1D("h_svLxySig_ee", "", 200, 0.0, 60);

   TH1D *h_electron_dR = new TH1D("h_electron_dR", "", 120, 0.0, 4.0);
   TH1D *h_kaon_ee_dR = new TH1D("h_kaon_ee_dR", "", 120, 0.0, 4.0);
   TH1D *h_jpsiPhiOpen_ee = new TH1D("h_jpsiPhiOpen_ee", "", 120, 0.0, 4.0);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if ((jentry%10000) == 0) std::cout<<"Procssing entry: "<<jentry<<std::endl;

      vector< pair<float,float> > jpsiPtList_ee;
      vector< pair<float,float> > phiPtList_ee;

      /***************************
       *
       * Electron Part
       *
       * ************************/

      for (int i=0; i<nEle; ++i) {

	 pair<float,float> jpsiPt_ee = make_pair(elePt_lead->at(i), elePt_sublead->at(i));
	 pair<float,float> phiPt_ee = make_pair(kaonEEPt_lead->at(i), kaonEEPt_sublead->at(i));

	 bool preCut = true;

	 if (eleCharge_lead->at(i) * eleCharge_sublead->at(i) < 0) continue;
	 if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) > 0) continue; // same charge
	 if (elePt_lead->at(i) < 0.4 || elePt_sublead->at(i) < 0.4 || kaonEEPt_lead->at(i) < 0.4 || kaonEEPt_sublead->at(i) < 0.4) continue;
         if (!eleConvVeto_lead->at(i) || !eleConvVeto_sublead->at(i)) continue;

	 if (preCut == true) {

	    if (find(jpsiPtList_ee.begin(), jpsiPtList_ee.end(), jpsiPt_ee) == jpsiPtList_ee.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_electron_pt_lead->Fill(elePt_lead->at(i));
	       h_electron_eta_lead->Fill(eleEta_lead->at(i));
	       h_electron_phi_lead->Fill(elePhi_lead->at(i));
	       h_electron_D0_lead->Fill(eleD0_lead->at(i));
	       h_electron_Dz_lead->Fill(eleDz_lead->at(i));

	       h_electron_pt_sublead->Fill(elePt_sublead->at(i));
	       h_electron_eta_sublead->Fill(eleEta_sublead->at(i));
	       h_electron_phi_sublead->Fill(elePhi_sublead->at(i));
	       h_electron_D0_sublead->Fill(eleD0_sublead->at(i));
	       h_electron_Dz_sublead->Fill(eleDz_sublead->at(i));

	       h_electron_dR->Fill(bsEEdRele->at(i));
	       h_q2_InvM->Fill(bsEEJpsiMass->at(i)*bsEEJpsiMass->at(i));

	       jpsiPtList_ee.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee.begin(), phiPtList_ee.end(), phiPt_ee) == phiPtList_ee.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_kaon_ee_pt_lead->Fill(kaonEEPt_lead->at(i));
	       h_kaon_ee_eta_lead->Fill(kaonEEEta_lead->at(i));
	       h_kaon_ee_phi_lead->Fill(kaonEEPhi_lead->at(i));
	       h_kaon_ee_Dz_lead->Fill(kaonEEDz_lead->at(i));
	       h_kaon_ee_D0_lead->Fill(kaonEED0_lead->at(i));
	       h_kaon_ee_D0Sig_lead->Fill(kaonEED0_lead->at(i)/kaonEED0Error_lead->at(i));
	       h_kaon_ee_normChi2_lead->Fill(kaonEETrkNormChi2_lead->at(i));

	       h_kaon_ee_pt_sublead->Fill(kaonEEPt_sublead->at(i));
	       h_kaon_ee_eta_sublead->Fill(kaonEEEta_sublead->at(i));
	       h_kaon_ee_phi_sublead->Fill(kaonEEPhi_sublead->at(i));
	       h_kaon_ee_Dz_sublead->Fill(kaonEEDz_sublead->at(i));
	       h_kaon_ee_D0_sublead->Fill(kaonEED0_sublead->at(i));
	       h_kaon_ee_D0Sig_sublead->Fill(kaonEED0_sublead->at(i)/kaonEED0Error_sublead->at(i));
	       h_kaon_ee_normChi2_sublead->Fill(kaonEETrkNormChi2_sublead->at(i));

	       h_kaon_ee_dR->Fill(bsEEdRkaon->at(i));
	       h_phiee_InvM->Fill(bsEEPhiMass->at(i));

	       phiPtList_ee.push_back(phiPt_ee);
	    }
	    //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	    //}
	    h_svProb_ee->Fill(eleSvProb->at(i));
	    h_svCosAngle_ee->Fill(eleSvCosAngle->at(i));
	    h_svCtxy_ee->Fill(eleSvCtxy->at(i));
	    h_svLxy_ee->Fill(eleSvLxy->at(i));
	    h_svLxySig_ee->Fill(eleSvLxy->at(i)/eleSvLxyError->at(i));

	    h_jpsiPhiOpen_ee->Fill(bsEEdRJpsiPhi->at(i));
	    h_bsee_InvM->Fill(bsEEBsMass->at(i));
	    

	 }

      }


   }

   file_out.cd();
   file_out.Write();
   file_out.Close();

}
