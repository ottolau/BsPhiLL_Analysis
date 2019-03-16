#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"

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

   // Histograms
   // Non-trigger Muon
   
   TH1D *h_muon_nontrg_pt_lead_preselected = new TH1D("h_muon_nontrg_pt_lead_preselected", "", 80, 0.0, 20.0);
   TH1D *h_muon_nontrg_pt_sublead_preselected = new TH1D("h_muon_nontrg_pt_sublead_preselected", "", 80, 0.0, 20.0);
   TH1D *h_muon_nontrg_eta_preselected = new TH1D("h_muon_nontrg_eta_preselected", "", 50, -3, 3);
   TH1D *h_muon_nontrg_phi_preselected = new TH1D("h_muon_nontrg_phi_preselected", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_muon_nontrg_dR_preselected = new TH1D("h_muon_nontrg_dR_preselected", "", 120, 0.0, 4.0);
   TH1D *h_muon_nontrg_D0_preselected = new TH1D("h_muon_nontrg_D0_preselected", "", 50, 0.0, 0.5);
   TH1D *h_muon_nontrg_Dz_preselected = new TH1D("h_muon_nontrg_Dz_preselected", "", 50, 0.0, 5.0);
   TH1D *h_muon_nontrg_SIP_preselected = new TH1D("h_muon_nontrg_SIP_preselected", "", 50, 0.0, 20.0);

   TH1D *h_kaon_uu_nontrg_pt_lead_preselected = new TH1D("h_kaon_uu_nontrg_pt_lead_preselected", "", 80, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_pt_sublead_preselected = new TH1D("h_kaon_uu_nontrg_pt_sublead_preselected", "", 80, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_eta_preselected = new TH1D("h_kaon_uu_nontrg_eta_preselected", "", 50, -3, 3);
   TH1D *h_kaon_uu_nontrg_phi_preselected = new TH1D("h_kaon_uu_nontrg_phi_preselected", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_uu_nontrg_dR_preselected = new TH1D("h_kaon_uu_nontrg_dR_preselected", "", 120, 0.0, 4.0);
   TH1D *h_kaon_uu_nontrg_D0_preselected = new TH1D("h_kaon_uu_nontrg_D0_preselected", "", 50, 0.0, 0.5);
   TH1D *h_kaon_uu_nontrg_Dz_preselected = new TH1D("h_kaon_uu_nontrg_Dz_preselected", "", 50, 0.0, 5.0);
   TH1D *h_kaon_uu_nontrg_D0Error_preselected = new TH1D("h_kaon_uu_nontrg_D0Error_preselected", "", 30, 0.0, 0.2);
   TH1D *h_kaon_uu_nontrg_DzError_preselected = new TH1D("h_kaon_uu_nontrg_DzError_preselected", "", 30, 0.0, 0.2);
   TH1D *h_kaon_uu_nontrg_D0Sig_preselected = new TH1D("h_kaon_uu_nontrg_D0Sig_preselected", "", 50, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_DzSig_preselected = new TH1D("h_kaon_uu_nontrg_DzSig_preselected", "", 100, 0.0, 100.0);
   TH1D *h_kaon_uu_nontrg_normChi2_preselected = new TH1D("h_kaon_uu_nontrg_normChi2_preselected", "", 100, 0.0, 5.0);

   TH1D *h_svCtxy_uu_nontrg_preselected = new TH1D("h_svCtxy_uu_nontrg_preselected", "", 100, 0.0, 200);
   TH1D *h_svChi2_uu_nontrg_preselected = new TH1D("h_svChi2_uu_nontrg_preselected", "", 100, 0.0, 5.0);
   TH1D *h_svProb_uu_nontrg_preselected = new TH1D("h_svProb_uu_nontrg_preselected", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_uu_nontrg_preselected = new TH1D("h_svCosAngle_uu_nontrg_preselected", "", 20, -1.0, 1.0);
   TH1D *h_jpsiPhiOpen_uu_nontrg_preselected = new TH1D("h_jpsiPhiOpen_uu_nontrg_preselected", "", 120, 0.0, 4.0);

   TH1D *h_muon_nontrg_pt_lead_selected = new TH1D("h_muon_nontrg_pt_lead_selected", "", 80, 0.0, 20.0);
   TH1D *h_muon_nontrg_pt_sublead_selected = new TH1D("h_muon_nontrg_pt_sublead_selected", "", 80, 0.0, 20.0);
   TH1D *h_muon_nontrg_eta_selected = new TH1D("h_muon_nontrg_eta_selected", "", 50, -3, 3);
   TH1D *h_muon_nontrg_phi_selected = new TH1D("h_muon_nontrg_phi_selected", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_muon_nontrg_dR_selected = new TH1D("h_muon_nontrg_dR_selected", "", 120, 0.0, 4.0);
   TH1D *h_muon_nontrg_D0_selected = new TH1D("h_muon_nontrg_D0_selected", "", 50, 0.0, 0.5);
   TH1D *h_muon_nontrg_Dz_selected = new TH1D("h_muon_nontrg_Dz_selected", "", 50, 0.0, 5.0);
   TH1D *h_muon_nontrg_SIP_selected = new TH1D("h_muon_nontrg_SIP_selected", "", 50, 0.0, 20.0);

   TH1D *h_kaon_uu_nontrg_pt_lead_selected = new TH1D("h_kaon_uu_nontrg_pt_lead_selected", "", 80, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_pt_sublead_selected = new TH1D("h_kaon_uu_nontrg_pt_sublead_selected", "", 80, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_eta_selected = new TH1D("h_kaon_uu_nontrg_eta_selected", "", 50, -3, 3);
   TH1D *h_kaon_uu_nontrg_phi_selected = new TH1D("h_kaon_uu_nontrg_phi_selected", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_kaon_uu_nontrg_dR_selected = new TH1D("h_kaon_uu_nontrg_dR_selected", "", 120, 0.0, 4.0);
   TH1D *h_kaon_uu_nontrg_D0_selected = new TH1D("h_kaon_uu_nontrg_D0_selected", "", 50, 0.0, 0.5);
   TH1D *h_kaon_uu_nontrg_Dz_selected = new TH1D("h_kaon_uu_nontrg_Dz_selected", "", 50, 0.0, 5.0);
   TH1D *h_kaon_uu_nontrg_D0Error_selected = new TH1D("h_kaon_uu_nontrg_D0Error_selected", "", 30, 0.0, 0.2);
   TH1D *h_kaon_uu_nontrg_DzError_selected = new TH1D("h_kaon_uu_nontrg_DzError_selected", "", 30, 0.0, 0.2);
   TH1D *h_kaon_uu_nontrg_D0Sig_selected = new TH1D("h_kaon_uu_nontrg_D0Sig_selected", "", 50, 0.0, 20.0);
   TH1D *h_kaon_uu_nontrg_DzSig_selected = new TH1D("h_kaon_uu_nontrg_DzSig_selected", "", 100, 0.0, 100.0);
   TH1D *h_kaon_uu_nontrg_normChi2_selected = new TH1D("h_kaon_uu_nontrg_normChi2_selected", "", 100, 0.0, 5.0);

   TH1D *h_svCtxy_uu_nontrg_selected = new TH1D("h_svCtxy_uu_nontrg_selected", "", 100, 0.0, 200);
   TH1D *h_svChi2_uu_nontrg_selected = new TH1D("h_svChi2_uu_nontrg_selected", "", 100, 0.0, 5.0);
   TH1D *h_svProb_uu_nontrg_selected = new TH1D("h_svProb_uu_nontrg_selected", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_uu_nontrg_selected = new TH1D("h_svCosAngle_uu_nontrg_selected", "", 20, -1.0, 1.0);
   TH1D *h_jpsiPhiOpen_uu_nontrg_selected = new TH1D("h_jpsiPhiOpen_uu_nontrg_selected", "", 120, 0.0, 4.0);

   TH1D *h_jpsiuu_nontrg_InvM_preselected = new TH1D("h_jpsiuu_nontrg_InvM_preselected", "", 140, 2.4, 3.8);
   TH1D *h_phiuu_nontrg_InvM_preselected = new TH1D("h_phiuu_nontrg_InvM_preselected", "", 480, 0.98, 1.1);
   TH1D *h_bsuu_nontrg_InvM_preselected = new TH1D("h_bsuu_nontrg_InvM_preselected", "", 400, 4.0, 8.0);
   TH1D *h_jpsiuu_nontrg_InvM_selected = new TH1D("h_jpsiuu_nontrg_InvM_selected", "", 140, 2.4, 3.8);
   TH1D *h_phiuu_nontrg_InvM_selected = new TH1D("h_phiuu_nontrg_InvM_selected", "", 480, 0.98, 1.1);
   TH1D *h_bsuu_nontrg_InvM_selected = new TH1D("h_bsuu_nontrg_InvM_selected", "", 400, 4.0, 8.0);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
 
      if ((jentry%100000) == 0) cout<<"Procssing entry: "<<jentry<<endl;
     
      vector< pair<float,float> > jpsiPtList_uu_nontrg_preselected;
      vector< pair<float,float> > phiPtList_uu_nontrg_preselected;
      vector< pair<float,float> > jpsiPtList_uu_nontrg_selected;
      vector< pair<float,float> > phiPtList_uu_nontrg_selected;

      /***************************
       *
       * Muon Part
       *
       * ************************/

      for (int i=0; i<nMu; ++i) {

	 pair<float,float> jpsiPt_uu = make_pair(muPt_lead->at(i), muPt_sublead->at(i));
	 pair<float,float> phiPt_uu = make_pair(kaonMMPt_lead->at(i), kaonMMPt_sublead->at(i));

	 // Non-trigger muon part
	 bool preCut_nontrg = true;

	 if (preCut_nontrg == true && (muFiredTrgs_lead->at(i) == false && muFiredTrgs_sublead->at(i) == false)) {

	    bool contam = false;
	    if (matchedTrg == false) {
	       for (int iTrg = 0; iTrg<nTrg; iTrg++){
		  if (deltaPhi(trgPhi->at(iTrg), muPhi_lead->at(i)) < 1.1 || deltaPhi(trgPhi->at(iTrg), muPhi_sublead->at(i)) < 1.1) {
		     contam = true;
		     break;
		  }
		  if (fabs(trgEta->at(iTrg)-muEta_lead->at(i)) < 0.6 || fabs(trgEta->at(iTrg)-muEta_sublead->at(i)) < 0.6) {
		     contam = true;
		     break;
		  }
	       }
	    }
	    
	    if (contam) continue;


	    if (find(jpsiPtList_uu_nontrg_preselected.begin(), jpsiPtList_uu_nontrg_preselected.end(), jpsiPt_uu) == jpsiPtList_uu_nontrg_preselected.end()) {
	       h_muon_nontrg_pt_lead_preselected->Fill(muPt_lead->at(i));
	       h_muon_nontrg_pt_sublead_preselected->Fill(muPt_sublead->at(i));
	       h_muon_nontrg_eta_preselected->Fill(muEta_lead->at(i));
	       h_muon_nontrg_phi_preselected->Fill(muPhi_lead->at(i));
	       h_muon_nontrg_D0_preselected->Fill(muD0_lead->at(i));
	       h_muon_nontrg_Dz_preselected->Fill(muDz_lead->at(i));
	       h_muon_nontrg_SIP_preselected->Fill(muSIP_lead->at(i));
	       h_muon_nontrg_eta_preselected->Fill(muEta_sublead->at(i));
	       h_muon_nontrg_phi_preselected->Fill(muPhi_sublead->at(i));
	       h_muon_nontrg_D0_preselected->Fill(muD0_sublead->at(i));
	       h_muon_nontrg_Dz_preselected->Fill(muDz_sublead->at(i));
	       h_muon_nontrg_SIP_preselected->Fill(muSIP_sublead->at(i));
	       h_muon_nontrg_dR_preselected->Fill(bsMMdRmu->at(i));
	       h_jpsiuu_nontrg_InvM_preselected->Fill(bsMMJpsiMass->at(i));

	       jpsiPtList_uu_nontrg_preselected.push_back(jpsiPt_uu);
	    }
	    if (find(phiPtList_uu_nontrg_preselected.begin(), phiPtList_uu_nontrg_preselected.end(), phiPt_uu) == phiPtList_uu_nontrg_preselected.end()) {
	       h_kaon_uu_nontrg_pt_lead_preselected->Fill(kaonMMPt_lead->at(i));
	       h_kaon_uu_nontrg_pt_sublead_preselected->Fill(kaonMMPt_sublead->at(i));
	       h_kaon_uu_nontrg_eta_preselected->Fill(kaonMMEta_lead->at(i));
	       h_kaon_uu_nontrg_phi_preselected->Fill(kaonMMPhi_lead->at(i));
	       h_kaon_uu_nontrg_D0_preselected->Fill(kaonMMD0_lead->at(i));
	       h_kaon_uu_nontrg_Dz_preselected->Fill(kaonMMDz_lead->at(i));
	       h_kaon_uu_nontrg_D0Error_preselected->Fill(kaonMMD0Error_lead->at(i));
	       h_kaon_uu_nontrg_DzError_preselected->Fill(kaonMMDzError_lead->at(i));
	       h_kaon_uu_nontrg_D0Sig_preselected->Fill(kaonMMD0_lead->at(i)/kaonMMD0Error_lead->at(i));
	       h_kaon_uu_nontrg_DzSig_preselected->Fill(kaonMMDz_lead->at(i)/kaonMMDzError_lead->at(i));
	       h_kaon_uu_nontrg_normChi2_preselected->Fill(kaonMMTrkNormChi2_lead->at(i));
	       h_kaon_uu_nontrg_eta_preselected->Fill(kaonMMEta_sublead->at(i));
	       h_kaon_uu_nontrg_phi_preselected->Fill(kaonMMPhi_sublead->at(i));
	       h_kaon_uu_nontrg_D0_preselected->Fill(kaonMMD0_sublead->at(i));
	       h_kaon_uu_nontrg_Dz_preselected->Fill(kaonMMDz_sublead->at(i));
	       h_kaon_uu_nontrg_D0Error_preselected->Fill(kaonMMD0Error_sublead->at(i));
	       h_kaon_uu_nontrg_DzError_preselected->Fill(kaonMMDzError_sublead->at(i));
	       h_kaon_uu_nontrg_D0Sig_preselected->Fill(kaonMMD0_sublead->at(i)/kaonMMD0Error_sublead->at(i));
	       h_kaon_uu_nontrg_DzSig_preselected->Fill(kaonMMDz_sublead->at(i)/kaonMMDzError_sublead->at(i));
	       h_kaon_uu_nontrg_normChi2_preselected->Fill(kaonMMTrkNormChi2_sublead->at(i));
	       h_kaon_uu_nontrg_dR_preselected->Fill(bsMMdRkaon->at(i));
	       h_phiuu_nontrg_InvM_preselected->Fill(bsMMPhiMass->at(i));

	       phiPtList_uu_nontrg_preselected.push_back(phiPt_uu);
	    }
	    h_svCtxy_uu_nontrg_preselected->Fill(muSvCtxy->at(i));
	    h_svChi2_uu_nontrg_preselected->Fill(muSvChi2->at(i)/muSvNDOF->at(i));
	    h_svProb_uu_nontrg_preselected->Fill(TMath::Prob(muSvChi2->at(i), muSvNDOF->at(i)));
	    h_svCosAngle_uu_nontrg_preselected->Fill(muSvCosAngle->at(i));
	    h_jpsiPhiOpen_uu_nontrg_preselected->Fill(bsMMdRJpsiPhi->at(i));

	    if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup) {
	       h_bsuu_nontrg_InvM_preselected->Fill(bsMMBsMass->at(i));

	    }

	    bool cutBasedSelected = true;
	    if (!(muPt_lead->at(i) > 2 && muPt_sublead->at(i) > 2)) cutBasedSelected = false;
	    if (!(kaonMMPt_lead->at(i) > 0.7 && kaonMMPt_sublead->at(i) > 0.7)) cutBasedSelected = false;


	    if (cutBasedSelected) {

	       if (find(jpsiPtList_uu_nontrg_selected.begin(), jpsiPtList_uu_nontrg_selected.end(), jpsiPt_uu) == jpsiPtList_uu_nontrg_selected.end()) {
		  h_muon_nontrg_pt_lead_selected->Fill(muPt_lead->at(i));
		  h_muon_nontrg_pt_sublead_selected->Fill(muPt_sublead->at(i));
		  h_muon_nontrg_eta_selected->Fill(muEta_lead->at(i));
		  h_muon_nontrg_phi_selected->Fill(muPhi_lead->at(i));
		  h_muon_nontrg_D0_selected->Fill(muD0_lead->at(i));
		  h_muon_nontrg_Dz_selected->Fill(muDz_lead->at(i));
		  h_muon_nontrg_SIP_selected->Fill(muSIP_lead->at(i));
		  h_muon_nontrg_eta_selected->Fill(muEta_sublead->at(i));
		  h_muon_nontrg_phi_selected->Fill(muPhi_sublead->at(i));
		  h_muon_nontrg_D0_selected->Fill(muD0_sublead->at(i));
		  h_muon_nontrg_Dz_selected->Fill(muDz_sublead->at(i));
		  h_muon_nontrg_SIP_selected->Fill(muSIP_sublead->at(i));
		  h_muon_nontrg_dR_selected->Fill(bsMMdRmu->at(i));
		  h_jpsiuu_nontrg_InvM_selected->Fill(bsMMJpsiMass->at(i));

		  jpsiPtList_uu_nontrg_selected.push_back(jpsiPt_uu);
	       }
	       if (find(phiPtList_uu_nontrg_selected.begin(), phiPtList_uu_nontrg_selected.end(), phiPt_uu) == phiPtList_uu_nontrg_selected.end()) {
		  h_kaon_uu_nontrg_pt_lead_selected->Fill(kaonMMPt_lead->at(i));
		  h_kaon_uu_nontrg_pt_sublead_selected->Fill(kaonMMPt_sublead->at(i));
		  h_kaon_uu_nontrg_eta_selected->Fill(kaonMMEta_lead->at(i));
		  h_kaon_uu_nontrg_phi_selected->Fill(kaonMMPhi_lead->at(i));
		  h_kaon_uu_nontrg_D0_selected->Fill(kaonMMD0_lead->at(i));
		  h_kaon_uu_nontrg_Dz_selected->Fill(kaonMMDz_lead->at(i));
		  h_kaon_uu_nontrg_D0Error_selected->Fill(kaonMMD0Error_lead->at(i));
		  h_kaon_uu_nontrg_DzError_selected->Fill(kaonMMDzError_lead->at(i));
		  h_kaon_uu_nontrg_D0Sig_selected->Fill(kaonMMD0_lead->at(i)/kaonMMD0Error_lead->at(i));
		  h_kaon_uu_nontrg_DzSig_selected->Fill(kaonMMDz_lead->at(i)/kaonMMDzError_lead->at(i));
		  h_kaon_uu_nontrg_normChi2_selected->Fill(kaonMMTrkNormChi2_lead->at(i));
		  h_kaon_uu_nontrg_eta_selected->Fill(kaonMMEta_sublead->at(i));
		  h_kaon_uu_nontrg_phi_selected->Fill(kaonMMPhi_sublead->at(i));
		  h_kaon_uu_nontrg_D0_selected->Fill(kaonMMD0_sublead->at(i));
		  h_kaon_uu_nontrg_Dz_selected->Fill(kaonMMDz_sublead->at(i));
		  h_kaon_uu_nontrg_D0Error_selected->Fill(kaonMMD0Error_sublead->at(i));
		  h_kaon_uu_nontrg_DzError_selected->Fill(kaonMMDzError_sublead->at(i));
		  h_kaon_uu_nontrg_D0Sig_selected->Fill(kaonMMD0_sublead->at(i)/kaonMMD0Error_sublead->at(i));
		  h_kaon_uu_nontrg_DzSig_selected->Fill(kaonMMDz_sublead->at(i)/kaonMMDzError_sublead->at(i));
		  h_kaon_uu_nontrg_normChi2_selected->Fill(kaonMMTrkNormChi2_sublead->at(i));
		  h_kaon_uu_nontrg_dR_selected->Fill(bsMMdRkaon->at(i));
		  h_phiuu_nontrg_InvM_selected->Fill(bsMMPhiMass->at(i));

		  phiPtList_uu_nontrg_selected.push_back(phiPt_uu);
	       }
	       h_svCtxy_uu_nontrg_selected->Fill(muSvCtxy->at(i));
	       h_svChi2_uu_nontrg_selected->Fill(muSvChi2->at(i)/muSvNDOF->at(i));
	       h_svProb_uu_nontrg_selected->Fill(TMath::Prob(muSvChi2->at(i), muSvNDOF->at(i)));
	       h_svCosAngle_uu_nontrg_selected->Fill(muSvCosAngle->at(i));
	       h_jpsiPhiOpen_uu_nontrg_selected->Fill(bsMMdRJpsiPhi->at(i));

	       if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup) {

		  h_bsuu_nontrg_InvM_selected->Fill(bsMMBsMass->at(i));

	       }

	    }

	 }

      }

   }

   file_out.cd();
   file_out.Write();
   file_out.Close();


}
