#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry > maxevents && maxevents != -1) break;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      vector< pair<float,float> > jpsiPtList_ee;
      vector< pair<float,float> > phiPtList_ee;
      vector< pair<float,float> > jpsiPtList_uu;
      vector< pair<float,float> > phiPtList_uu;
      vector< pair<float,float> > jpsiPtList_uu_nontrg;
      vector< pair<float,float> > phiPtList_uu_nontrg;

      /***************************
       *
       * Electron Part
       *
       * ************************/

      for (int i=0; i<nEle; ++i) {

	 pair<float,float> jpsiPt_ee = make_pair(elePt_lead->at(i), elePt_sublead->at(i));
	 pair<float,float> phiPt_ee = make_pair(kaonEEPt_lead->at(i), kaonEEPt_sublead->at(i));

	 bool preCut = false;

	 if (preCut == true) {

	    if (find(jpsiPtList_ee.begin(), jpsiPtList_ee.end(), jpsiPt_ee) == jpsiPtList_ee.end()) {
	       if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {

	       }

	       jpsiPtList_ee.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee.begin(), phiPtList_ee.end(), phiPt_ee) == phiPtList_ee.end()) {
	       if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {

	       }

	       phiPtList_ee.push_back(phiPt_ee);
	    }
	    if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {

	    }
	 }

      }

      /***************************
       *
       * Muon Part
       *
       * ************************/

      for (int i=0; i<nMu; ++i) {

	 pair<float,float> jpsiPt_uu = make_pair(muPt_lead->at(i), muPt_sublead->at(i));
	 pair<float,float> phiPt_uu = make_pair(kaonMMPt_lead->at(i), kaonMMPt_sublead->at(i));


	 // Trigger muon part
	 bool preCut_trgmu = false;

	 if (preCut_trgmu == true && (muFiredTrgs_lead->at(i) == true || muFiredTrgs_sublead->at(i) == true)) {

	    if (find(jpsiPtList_uu.begin(), jpsiPtList_uu.end(), jpsiPt_uu) == jpsiPtList_uu.end()) {
	       if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	       }

	       jpsiPtList_uu.push_back(jpsiPt_uu);
	    }
	    if (find(phiPtList_uu.begin(), phiPtList_uu.end(), phiPt_uu) == phiPtList_uu.end()) {
	       if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	       }

	       phiPtList_uu.push_back(phiPt_uu);
	    }
	    if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	    }
	 }

      


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


	    if (find(jpsiPtList_uu_nontrg.begin(), jpsiPtList_uu_nontrg.end(), jpsiPt_uu) == jpsiPtList_uu_nontrg.end()) {
	       if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	       }

	       jpsiPtList_uu_nontrg.push_back(jpsiPt_uu);
	    }
	    if (find(phiPtList_uu_nontrg.begin(), phiPtList_uu_nontrg.end(), phiPt_uu) == phiPtList_uu_nontrg.end()) {
	       if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	       }

	       phiPtList_uu_nontrg.push_back(phiPt_uu);
	    }
	    if (bsMMJpsiMass->at(i) > jpsilow && bsMMJpsiMass->at(i) < jpsiup && bsMMPhiMass->at(i) > philow && bsMMPhiMass->at(i) < phiup && bsMMBsMass->at(i) > bslow && bsMMBsMass->at(i) < bsup) {

	    }
	 }

      }


   }
}
