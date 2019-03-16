#define BsPhiLLTupleTree_cxx
#include "BsPhiLLTupleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"

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

   TMVA::Tools::Instance();
   TMVA::PyMethodBase::PyInitialize();

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   reader->AddVariable("elePtLead", &elePtLead);
   reader->AddVariable("elePtSublead", &elePtSublead);
   reader->AddVariable("kaonPtLead", &kaonPtLead);
   reader->AddVariable("kaonPtSublead", &kaonPtSublead);
   reader->AddVariable("jpsiPt", &jpsiPt);
   reader->AddVariable("phiPt", &phiPt);
   reader->AddVariable("bsPt", &bsPt);
   reader->AddVariable("eledR", &eledR);
   reader->AddVariable("kaondR", &kaondR);
   reader->AddVariable("jpsiPhidR", &jpsiPhidR);
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

   //TString TMVAmethod = "PyKeras method";
   TString TMVAmethod = "PyGTB method";

   reader->BookMVA("PyKeras method", "/uscms/home/klau/nobackup/BtoKll/BParking2018/BsPhiLL_Integrated/BsPhiLL_Analysis/MVATraining/dataset/weights/TMVAClassification_BsPhiEE_PyKeras.weights.xml");
   reader->BookMVA("PyGTB method", "/uscms/home/klau/nobackup/BtoKll/BParking2018/BsPhiLL_Integrated/BsPhiLL_Analysis/MVATraining/dataset/weights/TMVAClassification_BsPhiEE_PyGTB.weights.xml");


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

   TH1D *h_electron_pt_sublead_pass = new TH1D("h_electron_pt_sublead_pass", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_sublead_pass = new TH1D("h_electron_eta_sublead_pass", "", 50, -3, 3);
   TH1D *h_electron_phi_sublead_pass = new TH1D("h_electron_phi_sublead_pass", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_sublead_pass = new TH1D("h_electron_D0_sublead_pass", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_sublead_pass = new TH1D("h_electron_Dz_sublead_pass", "", 100, -1.0, 1.0);

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

   TH1D *h_svProb_ee_pass = new TH1D("h_svProb_ee_pass", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_ee_pass = new TH1D("h_svCosAngle_ee_pass", "", 20, -1.0, 1.0);
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

   TH1D *h_electron_pt_sublead_veto = new TH1D("h_electron_pt_sublead_veto", "", 80, 0.0, 20.0);
   TH1D *h_electron_eta_sublead_veto = new TH1D("h_electron_eta_sublead_veto", "", 50, -3, 3);
   TH1D *h_electron_phi_sublead_veto = new TH1D("h_electron_phi_sublead_veto", "", 50, -1.0*TMath::Pi(), TMath::Pi());
   TH1D *h_electron_D0_sublead_veto = new TH1D("h_electron_D0_sublead_veto", "", 100, -0.5, 0.5);
   TH1D *h_electron_Dz_sublead_veto = new TH1D("h_electron_Dz_sublead_veto", "", 100, -1.0, 1.0);

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

   TH1D *h_svProb_ee_veto = new TH1D("h_svProb_ee_veto", "", 100, 0.0, 1.0);
   TH1D *h_svCosAngle_ee_veto = new TH1D("h_svCosAngle_ee_veto", "", 20, -1.0, 1.0);
   TH1D *h_svCtxy_ee_veto = new TH1D("h_svCtxy_ee_veto", "", 100, -1.0, 1.0);
   TH1D *h_svLxy_ee_veto = new TH1D("h_svLxy_ee_veto", "", 100, 0.0, 1);
   TH1D *h_svLxySig_ee_veto = new TH1D("h_svLxySig_ee_veto", "", 200, 0.0, 60);

   TH1D *h_electron_dR_veto = new TH1D("h_electron_dR_veto", "", 120, 0.0, 4.0);
   TH1D *h_kaon_ee_dR_veto = new TH1D("h_kaon_ee_dR_veto", "", 120, 0.0, 4.0);
   TH1D *h_jpsiPhiOpen_ee_veto = new TH1D("h_jpsiPhiOpen_ee_veto", "", 120, 0.0, 4.0);

   TH1D *h_mvaValue = new TH1D("h_mvaValue", "", 100, 0.0, 1.0);


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

      for (int i=0; i<nEle; ++i) {

	 bool preCut = true;

	 if (eleCharge_lead->at(i) * eleCharge_sublead->at(i) < 0) continue;
	 if (kaonEECharge_lead->at(i) * kaonEECharge_sublead->at(i) < 0) continue; // same charge
	 if (elePt_lead->at(i) < 0.4 || elePt_sublead->at(i) < 0.4 || kaonEEPt_lead->at(i) < 0.4 || kaonEEPt_sublead->at(i) < 0.4) continue;
         if (!eleConvVeto_lead->at(i) || !eleConvVeto_sublead->at(i)) continue;

	 pair<float,float> jpsiPt_ee = make_pair(elePt_lead->at(i), elePt_sublead->at(i));
	 pair<float,float> phiPt_ee = make_pair(kaonEEPt_lead->at(i), kaonEEPt_sublead->at(i));

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

	 double mvaValue = reader->EvaluateMVA(TMVAmethod);

	 h_mvaValue->Fill(mvaValue);

	 if (mvaValue > mvaCut) {

	    if (find(jpsiPtList_ee_pass.begin(), jpsiPtList_ee_pass.end(), jpsiPt_ee) == jpsiPtList_ee_pass.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_electron_pt_lead_pass->Fill(elePt_lead->at(i));
	       h_electron_eta_lead_pass->Fill(eleEta_lead->at(i));
	       h_electron_phi_lead_pass->Fill(elePhi_lead->at(i));
	       h_electron_D0_lead_pass->Fill(eleD0_lead->at(i));
	       h_electron_Dz_lead_pass->Fill(eleDz_lead->at(i));

	       h_electron_pt_sublead_pass->Fill(elePt_sublead->at(i));
	       h_electron_eta_sublead_pass->Fill(eleEta_sublead->at(i));
	       h_electron_phi_sublead_pass->Fill(elePhi_sublead->at(i));
	       h_electron_D0_sublead_pass->Fill(eleD0_sublead->at(i));
	       h_electron_Dz_sublead_pass->Fill(eleDz_sublead->at(i));

	       h_electron_dR_pass->Fill(bsEEdRele->at(i));
	       h_q2_InvM_pass->Fill(bsEEJpsiMass->at(i)*bsEEJpsiMass->at(i));
	       h_jpsi_InvM_pass->Fill(bsEEJpsiMass->at(i));

	       jpsiPtList_ee_pass.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee_pass.begin(), phiPtList_ee_pass.end(), phiPt_ee) == phiPtList_ee_pass.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_kaon_ee_pt_lead_pass->Fill(kaonEEPt_lead->at(i));
	       h_kaon_ee_eta_lead_pass->Fill(kaonEEEta_lead->at(i));
	       h_kaon_ee_phi_lead_pass->Fill(kaonEEPhi_lead->at(i));
	       h_kaon_ee_Dz_lead_pass->Fill(kaonEEDz_lead->at(i));
	       h_kaon_ee_D0_lead_pass->Fill(kaonEED0_lead->at(i));
	       h_kaon_ee_D0Sig_lead_pass->Fill(kaonEED0_lead->at(i)/kaonEED0Error_lead->at(i));
	       h_kaon_ee_normChi2_lead_pass->Fill(kaonEETrkNormChi2_lead->at(i));

	       h_kaon_ee_pt_sublead_pass->Fill(kaonEEPt_sublead->at(i));
	       h_kaon_ee_eta_sublead_pass->Fill(kaonEEEta_sublead->at(i));
	       h_kaon_ee_phi_sublead_pass->Fill(kaonEEPhi_sublead->at(i));
	       h_kaon_ee_Dz_sublead_pass->Fill(kaonEEDz_sublead->at(i));
	       h_kaon_ee_D0_sublead_pass->Fill(kaonEED0_sublead->at(i));
	       h_kaon_ee_D0Sig_sublead_pass->Fill(kaonEED0_sublead->at(i)/kaonEED0Error_sublead->at(i));
	       h_kaon_ee_normChi2_sublead_pass->Fill(kaonEETrkNormChi2_sublead->at(i));

	       h_kaon_ee_dR_pass->Fill(bsEEdRkaon->at(i));
	       h_phiee_InvM_pass->Fill(bsEEPhiMass->at(i));

	       phiPtList_ee_pass.push_back(phiPt_ee);
	    }
	    //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	    //}
	    h_svProb_ee_pass->Fill(eleSvProb->at(i));
	    h_svCosAngle_ee_pass->Fill(eleSvCosAngle->at(i));
	    h_svCtxy_ee_pass->Fill(eleSvCtxy->at(i));
	    h_svLxy_ee_pass->Fill(eleSvLxy->at(i));
	    h_svLxySig_ee_pass->Fill(eleSvLxy->at(i)/eleSvLxyError->at(i));

	    h_jpsiPhiOpen_ee_pass->Fill(bsEEdRJpsiPhi->at(i));

	    if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup) {
	       h_bsee_InvM_pass->Fill(bsEEBsMass->at(i));
	    }

	 } else {

	    if (find(jpsiPtList_ee_veto.begin(), jpsiPtList_ee_veto.end(), jpsiPt_ee) == jpsiPtList_ee_veto.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_electron_pt_lead_veto->Fill(elePt_lead->at(i));
	       h_electron_eta_lead_veto->Fill(eleEta_lead->at(i));
	       h_electron_phi_lead_veto->Fill(elePhi_lead->at(i));
	       h_electron_D0_lead_veto->Fill(eleD0_lead->at(i));
	       h_electron_Dz_lead_veto->Fill(eleDz_lead->at(i));

	       h_electron_pt_sublead_veto->Fill(elePt_sublead->at(i));
	       h_electron_eta_sublead_veto->Fill(eleEta_sublead->at(i));
	       h_electron_phi_sublead_veto->Fill(elePhi_sublead->at(i));
	       h_electron_D0_sublead_veto->Fill(eleD0_sublead->at(i));
	       h_electron_Dz_sublead_veto->Fill(eleDz_sublead->at(i));

	       h_electron_dR_veto->Fill(bsEEdRele->at(i));
	       h_q2_InvM_veto->Fill(bsEEJpsiMass->at(i)*bsEEJpsiMass->at(i));
	       h_jpsi_InvM_veto->Fill(bsEEJpsiMass->at(i));

	       jpsiPtList_ee_veto.push_back(jpsiPt_ee);
	    }
	    if (find(phiPtList_ee_veto.begin(), phiPtList_ee_veto.end(), phiPt_ee) == phiPtList_ee_veto.end()) {
	       //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	       //}
	       h_kaon_ee_pt_lead_veto->Fill(kaonEEPt_lead->at(i));
	       h_kaon_ee_eta_lead_veto->Fill(kaonEEEta_lead->at(i));
	       h_kaon_ee_phi_lead_veto->Fill(kaonEEPhi_lead->at(i));
	       h_kaon_ee_Dz_lead_veto->Fill(kaonEEDz_lead->at(i));
	       h_kaon_ee_D0_lead_veto->Fill(kaonEED0_lead->at(i));
	       h_kaon_ee_D0Sig_lead_veto->Fill(kaonEED0_lead->at(i)/kaonEED0Error_lead->at(i));
	       h_kaon_ee_normChi2_lead_veto->Fill(kaonEETrkNormChi2_lead->at(i));

	       h_kaon_ee_pt_sublead_veto->Fill(kaonEEPt_sublead->at(i));
	       h_kaon_ee_eta_sublead_veto->Fill(kaonEEEta_sublead->at(i));
	       h_kaon_ee_phi_sublead_veto->Fill(kaonEEPhi_sublead->at(i));
	       h_kaon_ee_Dz_sublead_veto->Fill(kaonEEDz_sublead->at(i));
	       h_kaon_ee_D0_sublead_veto->Fill(kaonEED0_sublead->at(i));
	       h_kaon_ee_D0Sig_sublead_veto->Fill(kaonEED0_sublead->at(i)/kaonEED0Error_sublead->at(i));
	       h_kaon_ee_normChi2_sublead_veto->Fill(kaonEETrkNormChi2_sublead->at(i));

	       h_kaon_ee_dR_veto->Fill(bsEEdRkaon->at(i));
	       h_phiee_InvM_veto->Fill(bsEEPhiMass->at(i));

	       phiPtList_ee_veto.push_back(phiPt_ee);
	    }
	    //if (bsEEJpsiMass->at(i) > jpsilow && bsEEJpsiMass->at(i) < jpsiup && bsEEPhiMass->at(i) > philow && bsEEPhiMass->at(i) < phiup && bsEEBsMass->at(i) > bslow && bsEEBsMass->at(i) < bsup) {
	    //}
	    h_svProb_ee_veto->Fill(eleSvProb->at(i));
	    h_svCosAngle_ee_veto->Fill(eleSvCosAngle->at(i));
	    h_svCtxy_ee_veto->Fill(eleSvCtxy->at(i));
	    h_svLxy_ee_veto->Fill(eleSvLxy->at(i));
	    h_svLxySig_ee_veto->Fill(eleSvLxy->at(i)/eleSvLxyError->at(i));

	    h_jpsiPhiOpen_ee_veto->Fill(bsEEdRJpsiPhi->at(i));
	    h_bsee_InvM_veto->Fill(bsEEBsMass->at(i));
	    

	 }

      }


   }

   file_out.cd();
   file_out.Write();
   file_out.Close();

}
