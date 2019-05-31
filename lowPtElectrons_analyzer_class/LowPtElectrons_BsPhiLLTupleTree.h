//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  3 11:03:29 2019 by ROOT version 6.14/04
// from TTree EventTree/Event data (tag V08_00_26_03)
// found on file: ggtree_data.root
//////////////////////////////////////////////////////////

#ifndef LowPtElectrons_BsPhiLLTupleTree_h
#define LowPtElectrons_BsPhiLLTupleTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TH1.h"
#include "TMath.h"
#include "vector"
#include "math.h"

class LowPtElectrons_BsPhiLLTupleTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Int_t           nTrksPV;
   Bool_t          isPVGood;
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTPhoRejectedByPS;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   vector<int>     *phoPrescale;
   Int_t           nTrg;
   vector<float>   *trgPt;
   vector<float>   *trgEta;
   vector<float>   *trgPhi;
   Int_t           hltMu9IP6;
   vector<string>  *trgPath;
   Int_t           nLowPt;
   vector<int>     *lowPtCharge_lead;
   vector<float>   *lowPtD0_lead;
   vector<float>   *lowPtDz_lead;
   vector<float>   *lowPtD0Error_lead;
   vector<float>   *lowPtDzError_lead;
   vector<float>   *lowPtPt_lead;
   vector<float>   *lowPtEta_lead;
   vector<float>   *lowPtPhi_lead;
   vector<float>   *lowPtMVABWP_lead;
   vector<float>   *lowPtMVAUnBWP_lead;
   vector<int>     *lowPtCharge_sublead;
   vector<float>   *lowPtD0_sublead;
   vector<float>   *lowPtDz_sublead;
   vector<float>   *lowPtD0Error_sublead;
   vector<float>   *lowPtDzError_sublead;
   vector<float>   *lowPtPt_sublead;
   vector<float>   *lowPtEta_sublead;
   vector<float>   *lowPtPhi_sublead;
   vector<float>   *lowPtMVABWP_sublead;
   vector<float>   *lowPtMVAUnBWP_sublead;
   vector<float>   *lowPtSvChi2;
   vector<float>   *lowPtSvNDOF;
   vector<float>   *lowPtSvProb;
   vector<float>   *lowPtSvX;
   vector<float>   *lowPtSvY;
   vector<float>   *lowPtSvZ;
   vector<float>   *lowPtSvXError;
   vector<float>   *lowPtSvYError;
   vector<float>   *lowPtSvZError;
   vector<float>   *lowPtSvMass;
   vector<float>   *lowPtSvCtxy;
   vector<float>   *lowPtSvCosAngle;
   vector<float>   *lowPtSvLxy;
   vector<float>   *lowPtSvLxyError;
   vector<int>     *kaonLowPtCharge_lead;
   vector<float>   *kaonLowPtD0_lead;
   vector<float>   *kaonLowPtDz_lead;
   vector<float>   *kaonLowPtD0Error_lead;
   vector<float>   *kaonLowPtDzError_lead;
   vector<float>   *kaonLowPtPt_lead;
   vector<float>   *kaonLowPtEta_lead;
   vector<float>   *kaonLowPtPhi_lead;
   vector<float>   *kaonLowPtVx_lead;
   vector<float>   *kaonLowPtVy_lead;
   vector<float>   *kaonLowPtVz_lead;
   vector<float>   *kaonLowPtTrkChi2_lead;
   vector<float>   *kaonLowPtTrkNDOF_lead;
   vector<float>   *kaonLowPtTrkNormChi2_lead;
   vector<int>     *kaonLowPtCharge_sublead;
   vector<float>   *kaonLowPtD0_sublead;
   vector<float>   *kaonLowPtDz_sublead;
   vector<float>   *kaonLowPtD0Error_sublead;
   vector<float>   *kaonLowPtDzError_sublead;
   vector<float>   *kaonLowPtPt_sublead;
   vector<float>   *kaonLowPtEta_sublead;
   vector<float>   *kaonLowPtPhi_sublead;
   vector<float>   *kaonLowPtVx_sublead;
   vector<float>   *kaonLowPtVy_sublead;
   vector<float>   *kaonLowPtVz_sublead;
   vector<float>   *kaonLowPtTrkChi2_sublead;
   vector<float>   *kaonLowPtTrkNDOF_sublead;
   vector<float>   *kaonLowPtTrkNormChi2_sublead;
   vector<float>   *bsLowPtdRele;
   vector<float>   *bsLowPtdRkaon;
   vector<float>   *bsLowPtdRJpsiPhi;
   vector<float>   *bsLowPtJpsiMass;
   vector<float>   *bsLowPtPhiMass;
   vector<float>   *bsLowPtBsMass;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTPhoRejectedByPS;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_phoPrescale;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_trgPt;   //!
   TBranch        *b_trgEta;   //!
   TBranch        *b_trgPhi;   //!
   TBranch        *b_hltMu9IP6;   //!
   TBranch        *b_trgPath;   //!
   TBranch        *b_nLowPt;   //!
   TBranch        *b_lowPtCharge_lead;   //!
   TBranch        *b_lowPtD0_lead;   //!
   TBranch        *b_lowPtDz_lead;   //!
   TBranch        *b_lowPtD0Error_lead;   //!
   TBranch        *b_lowPtDzError_lead;   //!
   TBranch        *b_lowPtPt_lead;   //!
   TBranch        *b_lowPtEta_lead;   //!
   TBranch        *b_lowPtPhi_lead;   //!
   TBranch        *b_lowPtMVABWP_lead;   //!
   TBranch        *b_lowPtMVAUnBWP_lead;   //!
   TBranch        *b_lowPtCharge_sublead;   //!
   TBranch        *b_lowPtD0_sublead;   //!
   TBranch        *b_lowPtDz_sublead;   //!
   TBranch        *b_lowPtD0Error_sublead;   //!
   TBranch        *b_lowPtDzError_sublead;   //!
   TBranch        *b_lowPtPt_sublead;   //!
   TBranch        *b_lowPtEta_sublead;   //!
   TBranch        *b_lowPtPhi_sublead;   //!
   TBranch        *b_lowPtMVABWP_sublead;   //!
   TBranch        *b_lowPtMVAUnBWP_sublead;   //!
   TBranch        *b_lowPtSvChi2;   //!
   TBranch        *b_lowPtSvNDOF;   //!
   TBranch        *b_lowPtSvProb;   //!
   TBranch        *b_lowPtSvX;   //!
   TBranch        *b_lowPtSvY;   //!
   TBranch        *b_lowPtSvZ;   //!
   TBranch        *b_lowPtSvXError;   //!
   TBranch        *b_lowPtSvYError;   //!
   TBranch        *b_lowPtSvZError;   //!
   TBranch        *b_lowPtSvMass;   //!
   TBranch        *b_lowPtSvCtxy;   //!
   TBranch        *b_lowPtSvCosAngle;   //!
   TBranch        *b_lowPtSvLxy;   //!
   TBranch        *b_lowPtSvLxyError;   //!
   TBranch        *b_kaonLowPtCharge_lead;   //!
   TBranch        *b_kaonLowPtD0_lead;   //!
   TBranch        *b_kaonLowPtDz_lead;   //!
   TBranch        *b_kaonLowPtD0Error_lead;   //!
   TBranch        *b_kaonLowPtDzError_lead;   //!
   TBranch        *b_kaonLowPtPt_lead;   //!
   TBranch        *b_kaonLowPtEta_lead;   //!
   TBranch        *b_kaonLowPtPhi_lead;   //!
   TBranch        *b_kaonLowPtVx_lead;   //!
   TBranch        *b_kaonLowPtVy_lead;   //!
   TBranch        *b_kaonLowPtVz_lead;   //!
   TBranch        *b_kaonLowPtTrkChi2_lead;   //!
   TBranch        *b_kaonLowPtTrkNDOF_lead;   //!
   TBranch        *b_kaonLowPtTrkNormChi2_lead;   //!
   TBranch        *b_kaonLowPtCharge_sublead;   //!
   TBranch        *b_kaonLowPtD0_sublead;   //!
   TBranch        *b_kaonLowPtDz_sublead;   //!
   TBranch        *b_kaonLowPtD0Error_sublead;   //!
   TBranch        *b_kaonLowPtDzError_sublead;   //!
   TBranch        *b_kaonLowPtPt_sublead;   //!
   TBranch        *b_kaonLowPtEta_sublead;   //!
   TBranch        *b_kaonLowPtPhi_sublead;   //!
   TBranch        *b_kaonLowPtVx_sublead;   //!
   TBranch        *b_kaonLowPtVy_sublead;   //!
   TBranch        *b_kaonLowPtVz_sublead;   //!
   TBranch        *b_kaonLowPtTrkChi2_sublead;   //!
   TBranch        *b_kaonLowPtTrkNDOF_sublead;   //!
   TBranch        *b_kaonLowPtTrkNormChi2_sublead;   //!
   TBranch        *b_bsLowPtdRele;   //!
   TBranch        *b_bsLowPtdRkaon;   //!
   TBranch        *b_bsLowPtdRJpsiPhi;   //!
   TBranch        *b_bsLowPtJpsiMass;   //!
   TBranch        *b_bsLowPtPhiMass;   //!
   TBranch        *b_bsLowPtBsMass;   //!

   Float_t elecM = 0.000510998;
   Float_t muonM = 0.1056583745;
   Float_t kaonM = 0.493677;
   Float_t jpsiM = 3.096916;
   Float_t phiM = 1.019461;
   Float_t bsM = 5.3663;
   Float_t jpsilow = 2.9;
   Float_t jpsiup = 3.2;
   Float_t philow = 1.01;
   Float_t phiup = 1.03;
   Float_t bslow = 5.0;
   Float_t bsup = 5.8;
   Float_t bssideminuslow = 5.1;
   Float_t bssideminusup = 5.3;
   Float_t bssidepluslow = 5.45;
   Float_t bssideplusup = 5.6;
   Float_t jpsisideminuslow = 2.4;
   Float_t jpsisideminusup = 2.8;
   Float_t jpsisidepluslow = 3.3;
   Float_t jpsisideplusup = 3.8;
   Float_t phisideminuslow = 1.0;
   Float_t phisideminusup = 1.01;
   Float_t phisidepluslow = 1.04;
   Float_t phisideplusup = 1.05;

   LowPtElectrons_BsPhiLLTupleTree(TTree *tree=0);
   virtual ~LowPtElectrons_BsPhiLLTupleTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outputfile, Int_t maxevents, Float_t mvaCut, bool oppCharge);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Float_t          deltaPhi(Float_t phi1, Float_t phi2);
};

#endif

#ifdef LowPtElectrons_BsPhiLLTupleTree_cxx
LowPtElectrons_BsPhiLLTupleTree::LowPtElectrons_BsPhiLLTupleTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ggtree_data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ggtree_data.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ggtree_data.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

LowPtElectrons_BsPhiLLTupleTree::~LowPtElectrons_BsPhiLLTupleTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LowPtElectrons_BsPhiLLTupleTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LowPtElectrons_BsPhiLLTupleTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void LowPtElectrons_BsPhiLLTupleTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   phoPrescale = 0;
   trgPt = 0;
   trgEta = 0;
   trgPhi = 0;
   trgPath = 0;
   lowPtCharge_lead = 0;
   lowPtD0_lead = 0;
   lowPtDz_lead = 0;
   lowPtD0Error_lead = 0;
   lowPtDzError_lead = 0;
   lowPtPt_lead = 0;
   lowPtEta_lead = 0;
   lowPtPhi_lead = 0;
   lowPtMVABWP_lead = 0;
   lowPtMVAUnBWP_lead = 0;
   lowPtCharge_sublead = 0;
   lowPtD0_sublead = 0;
   lowPtDz_sublead = 0;
   lowPtD0Error_sublead = 0;
   lowPtDzError_sublead = 0;
   lowPtPt_sublead = 0;
   lowPtEta_sublead = 0;
   lowPtPhi_sublead = 0;
   lowPtMVABWP_sublead = 0;
   lowPtMVAUnBWP_sublead = 0;
   lowPtSvChi2 = 0;
   lowPtSvNDOF = 0;
   lowPtSvProb = 0;
   lowPtSvX = 0;
   lowPtSvY = 0;
   lowPtSvZ = 0;
   lowPtSvXError = 0;
   lowPtSvYError = 0;
   lowPtSvZError = 0;
   lowPtSvMass = 0;
   lowPtSvCtxy = 0;
   lowPtSvCosAngle = 0;
   lowPtSvLxy = 0;
   lowPtSvLxyError = 0;
   kaonLowPtCharge_lead = 0;
   kaonLowPtD0_lead = 0;
   kaonLowPtDz_lead = 0;
   kaonLowPtD0Error_lead = 0;
   kaonLowPtDzError_lead = 0;
   kaonLowPtPt_lead = 0;
   kaonLowPtEta_lead = 0;
   kaonLowPtPhi_lead = 0;
   kaonLowPtVx_lead = 0;
   kaonLowPtVy_lead = 0;
   kaonLowPtVz_lead = 0;
   kaonLowPtTrkChi2_lead = 0;
   kaonLowPtTrkNDOF_lead = 0;
   kaonLowPtTrkNormChi2_lead = 0;
   kaonLowPtCharge_sublead = 0;
   kaonLowPtD0_sublead = 0;
   kaonLowPtDz_sublead = 0;
   kaonLowPtD0Error_sublead = 0;
   kaonLowPtDzError_sublead = 0;
   kaonLowPtPt_sublead = 0;
   kaonLowPtEta_sublead = 0;
   kaonLowPtPhi_sublead = 0;
   kaonLowPtVx_sublead = 0;
   kaonLowPtVy_sublead = 0;
   kaonLowPtVz_sublead = 0;
   kaonLowPtTrkChi2_sublead = 0;
   kaonLowPtTrkNDOF_sublead = 0;
   kaonLowPtTrkNormChi2_sublead = 0;
   bsLowPtdRele = 0;
   bsLowPtdRkaon = 0;
   bsLowPtdRJpsiPhi = 0;
   bsLowPtJpsiMass = 0;
   bsLowPtPhiMass = 0;
   bsLowPtBsMass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTPhoRejectedByPS", &HLTPhoRejectedByPS, &b_HLTPhoRejectedByPS);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("phoPrescale", &phoPrescale, &b_phoPrescale);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("trgPt", &trgPt, &b_trgPt);
   fChain->SetBranchAddress("trgEta", &trgEta, &b_trgEta);
   fChain->SetBranchAddress("trgPhi", &trgPhi, &b_trgPhi);
   fChain->SetBranchAddress("hltMu9IP6", &hltMu9IP6, &b_hltMu9IP6);
   fChain->SetBranchAddress("trgPath", &trgPath, &b_trgPath);
   fChain->SetBranchAddress("nLowPt", &nLowPt, &b_nLowPt);
   fChain->SetBranchAddress("lowPtCharge_lead", &lowPtCharge_lead, &b_lowPtCharge_lead);
   fChain->SetBranchAddress("lowPtD0_lead", &lowPtD0_lead, &b_lowPtD0_lead);
   fChain->SetBranchAddress("lowPtDz_lead", &lowPtDz_lead, &b_lowPtDz_lead);
   fChain->SetBranchAddress("lowPtD0Error_lead", &lowPtD0Error_lead, &b_lowPtD0Error_lead);
   fChain->SetBranchAddress("lowPtDzError_lead", &lowPtDzError_lead, &b_lowPtDzError_lead);
   fChain->SetBranchAddress("lowPtPt_lead", &lowPtPt_lead, &b_lowPtPt_lead);
   fChain->SetBranchAddress("lowPtEta_lead", &lowPtEta_lead, &b_lowPtEta_lead);
   fChain->SetBranchAddress("lowPtPhi_lead", &lowPtPhi_lead, &b_lowPtPhi_lead);
   fChain->SetBranchAddress("lowPtMVABWP_lead", &lowPtMVABWP_lead, &b_lowPtMVABWP_lead);
   fChain->SetBranchAddress("lowPtMVAUnBWP_lead", &lowPtMVAUnBWP_lead, &b_lowPtMVAUnBWP_lead);
   fChain->SetBranchAddress("lowPtCharge_sublead", &lowPtCharge_sublead, &b_lowPtCharge_sublead);
   fChain->SetBranchAddress("lowPtD0_sublead", &lowPtD0_sublead, &b_lowPtD0_sublead);
   fChain->SetBranchAddress("lowPtDz_sublead", &lowPtDz_sublead, &b_lowPtDz_sublead);
   fChain->SetBranchAddress("lowPtD0Error_sublead", &lowPtD0Error_sublead, &b_lowPtD0Error_sublead);
   fChain->SetBranchAddress("lowPtDzError_sublead", &lowPtDzError_sublead, &b_lowPtDzError_sublead);
   fChain->SetBranchAddress("lowPtPt_sublead", &lowPtPt_sublead, &b_lowPtPt_sublead);
   fChain->SetBranchAddress("lowPtEta_sublead", &lowPtEta_sublead, &b_lowPtEta_sublead);
   fChain->SetBranchAddress("lowPtPhi_sublead", &lowPtPhi_sublead, &b_lowPtPhi_sublead);
   fChain->SetBranchAddress("lowPtMVABWP_sublead", &lowPtMVABWP_sublead, &b_lowPtMVABWP_sublead);
   fChain->SetBranchAddress("lowPtMVAUnBWP_sublead", &lowPtMVAUnBWP_sublead, &b_lowPtMVAUnBWP_sublead);
   fChain->SetBranchAddress("lowPtSvChi2", &lowPtSvChi2, &b_lowPtSvChi2);
   fChain->SetBranchAddress("lowPtSvNDOF", &lowPtSvNDOF, &b_lowPtSvNDOF);
   fChain->SetBranchAddress("lowPtSvProb", &lowPtSvProb, &b_lowPtSvProb);
   fChain->SetBranchAddress("lowPtSvX", &lowPtSvX, &b_lowPtSvX);
   fChain->SetBranchAddress("lowPtSvY", &lowPtSvY, &b_lowPtSvY);
   fChain->SetBranchAddress("lowPtSvZ", &lowPtSvZ, &b_lowPtSvZ);
   fChain->SetBranchAddress("lowPtSvXError", &lowPtSvXError, &b_lowPtSvXError);
   fChain->SetBranchAddress("lowPtSvYError", &lowPtSvYError, &b_lowPtSvYError);
   fChain->SetBranchAddress("lowPtSvZError", &lowPtSvZError, &b_lowPtSvZError);
   fChain->SetBranchAddress("lowPtSvMass", &lowPtSvMass, &b_lowPtSvMass);
   fChain->SetBranchAddress("lowPtSvCtxy", &lowPtSvCtxy, &b_lowPtSvCtxy);
   fChain->SetBranchAddress("lowPtSvCosAngle", &lowPtSvCosAngle, &b_lowPtSvCosAngle);
   fChain->SetBranchAddress("lowPtSvLxy", &lowPtSvLxy, &b_lowPtSvLxy);
   fChain->SetBranchAddress("lowPtSvLxyError", &lowPtSvLxyError, &b_lowPtSvLxyError);
   fChain->SetBranchAddress("kaonLowPtCharge_lead", &kaonLowPtCharge_lead, &b_kaonLowPtCharge_lead);
   fChain->SetBranchAddress("kaonLowPtD0_lead", &kaonLowPtD0_lead, &b_kaonLowPtD0_lead);
   fChain->SetBranchAddress("kaonLowPtDz_lead", &kaonLowPtDz_lead, &b_kaonLowPtDz_lead);
   fChain->SetBranchAddress("kaonLowPtD0Error_lead", &kaonLowPtD0Error_lead, &b_kaonLowPtD0Error_lead);
   fChain->SetBranchAddress("kaonLowPtDzError_lead", &kaonLowPtDzError_lead, &b_kaonLowPtDzError_lead);
   fChain->SetBranchAddress("kaonLowPtPt_lead", &kaonLowPtPt_lead, &b_kaonLowPtPt_lead);
   fChain->SetBranchAddress("kaonLowPtEta_lead", &kaonLowPtEta_lead, &b_kaonLowPtEta_lead);
   fChain->SetBranchAddress("kaonLowPtPhi_lead", &kaonLowPtPhi_lead, &b_kaonLowPtPhi_lead);
   fChain->SetBranchAddress("kaonLowPtVx_lead", &kaonLowPtVx_lead, &b_kaonLowPtVx_lead);
   fChain->SetBranchAddress("kaonLowPtVy_lead", &kaonLowPtVy_lead, &b_kaonLowPtVy_lead);
   fChain->SetBranchAddress("kaonLowPtVz_lead", &kaonLowPtVz_lead, &b_kaonLowPtVz_lead);
   fChain->SetBranchAddress("kaonLowPtTrkChi2_lead", &kaonLowPtTrkChi2_lead, &b_kaonLowPtTrkChi2_lead);
   fChain->SetBranchAddress("kaonLowPtTrkNDOF_lead", &kaonLowPtTrkNDOF_lead, &b_kaonLowPtTrkNDOF_lead);
   fChain->SetBranchAddress("kaonLowPtTrkNormChi2_lead", &kaonLowPtTrkNormChi2_lead, &b_kaonLowPtTrkNormChi2_lead);
   fChain->SetBranchAddress("kaonLowPtCharge_sublead", &kaonLowPtCharge_sublead, &b_kaonLowPtCharge_sublead);
   fChain->SetBranchAddress("kaonLowPtD0_sublead", &kaonLowPtD0_sublead, &b_kaonLowPtD0_sublead);
   fChain->SetBranchAddress("kaonLowPtDz_sublead", &kaonLowPtDz_sublead, &b_kaonLowPtDz_sublead);
   fChain->SetBranchAddress("kaonLowPtD0Error_sublead", &kaonLowPtD0Error_sublead, &b_kaonLowPtD0Error_sublead);
   fChain->SetBranchAddress("kaonLowPtDzError_sublead", &kaonLowPtDzError_sublead, &b_kaonLowPtDzError_sublead);
   fChain->SetBranchAddress("kaonLowPtPt_sublead", &kaonLowPtPt_sublead, &b_kaonLowPtPt_sublead);
   fChain->SetBranchAddress("kaonLowPtEta_sublead", &kaonLowPtEta_sublead, &b_kaonLowPtEta_sublead);
   fChain->SetBranchAddress("kaonLowPtPhi_sublead", &kaonLowPtPhi_sublead, &b_kaonLowPtPhi_sublead);
   fChain->SetBranchAddress("kaonLowPtVx_sublead", &kaonLowPtVx_sublead, &b_kaonLowPtVx_sublead);
   fChain->SetBranchAddress("kaonLowPtVy_sublead", &kaonLowPtVy_sublead, &b_kaonLowPtVy_sublead);
   fChain->SetBranchAddress("kaonLowPtVz_sublead", &kaonLowPtVz_sublead, &b_kaonLowPtVz_sublead);
   fChain->SetBranchAddress("kaonLowPtTrkChi2_sublead", &kaonLowPtTrkChi2_sublead, &b_kaonLowPtTrkChi2_sublead);
   fChain->SetBranchAddress("kaonLowPtTrkNDOF_sublead", &kaonLowPtTrkNDOF_sublead, &b_kaonLowPtTrkNDOF_sublead);
   fChain->SetBranchAddress("kaonLowPtTrkNormChi2_sublead", &kaonLowPtTrkNormChi2_sublead, &b_kaonLowPtTrkNormChi2_sublead);
   fChain->SetBranchAddress("bsLowPtdRele", &bsLowPtdRele, &b_bsLowPtdRele);
   fChain->SetBranchAddress("bsLowPtdRkaon", &bsLowPtdRkaon, &b_bsLowPtdRkaon);
   fChain->SetBranchAddress("bsLowPtdRJpsiPhi", &bsLowPtdRJpsiPhi, &b_bsLowPtdRJpsiPhi);
   fChain->SetBranchAddress("bsLowPtJpsiMass", &bsLowPtJpsiMass, &b_bsLowPtJpsiMass);
   fChain->SetBranchAddress("bsLowPtPhiMass", &bsLowPtPhiMass, &b_bsLowPtPhiMass);
   fChain->SetBranchAddress("bsLowPtBsMass", &bsLowPtBsMass, &b_bsLowPtBsMass);
   Notify();
}

Bool_t LowPtElectrons_BsPhiLLTupleTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LowPtElectrons_BsPhiLLTupleTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LowPtElectrons_BsPhiLLTupleTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
Float_t LowPtElectrons_BsPhiLLTupleTree::deltaPhi(Float_t phi1, Float_t phi2)
{
   Float_t dphi = phi1 - phi2;
   while (dphi >= TMath::Pi()) dphi -= 2*TMath::Pi();
   while (dphi < -1*TMath::Pi()) dphi += 2*TMath::Pi();
   return dphi;
}
#endif // #ifdef LowPtElectrons_BsPhiLLTupleTree_cxx
