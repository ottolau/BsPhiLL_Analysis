//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 15 15:24:14 2018 by ROOT version 6.12/07
// from TTree EventTree/Event data (tag V08_00_26_03)
// found on file: ggtree_data.root
//////////////////////////////////////////////////////////

#ifndef BsPhiLLTupleTree_h
#define BsPhiLLTupleTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TH1.h"
#include "TMath.h"
#include "vector"
#include "math.h"

class BsPhiLLTupleTree {
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
   //Int_t           npfPho;
   //vector<float>   *pfphoEt;
   //vector<float>   *pfphoEta;
   //vector<float>   *pfphoPhi;
   Int_t           nEle;
   vector<int>     *eleCharge_lead;
   vector<int>     *eleChargeConsistent_lead;
   vector<float>   *eleEn_lead;
   vector<float>   *eleSCEn_lead;
   vector<float>   *eleEcalEn_lead;
   vector<float>   *eleESEnP1_lead;
   vector<float>   *eleESEnP2_lead;
   vector<float>   *eleD0_lead;
   vector<float>   *eleDz_lead;
   vector<float>   *eleD0Error_lead;
   vector<float>   *eleDzError_lead;
   vector<float>   *eleSIP_lead;
   vector<float>   *elePt_lead;
   vector<float>   *eleEta_lead;
   vector<float>   *elePhi_lead;
   vector<float>   *eleR9_lead;
   vector<float>   *eleCalibPt_lead;
   vector<float>   *eleCalibEn_lead;
   vector<float>   *eleSCEta_lead;
   vector<float>   *eleSCPhi_lead;
   vector<float>   *eleSCRawEn_lead;
   vector<float>   *eleSCEtaWidth_lead;
   vector<float>   *eleSCPhiWidth_lead;
   vector<float>   *eleHoverE_lead;
   vector<float>   *eleEoverP_lead;
   vector<float>   *eleEoverPout_lead;
   vector<float>   *eleEoverPInv_lead;
   vector<float>   *eleBrem_lead;
   vector<float>   *eledEtaAtVtx_lead;
   vector<float>   *eledPhiAtVtx_lead;
   vector<float>   *eledEtaAtCalo_lead;
   vector<float>   *eleSigmaIEtaIEtaFull5x5_lead;
   vector<float>   *eleSigmaIPhiIPhiFull5x5_lead;
   vector<int>     *eleConvVeto_lead;
   vector<int>     *eleMissHits_lead;
   vector<float>   *eleESEffSigmaRR_lead;
   vector<float>   *elePFChIso_lead;
   vector<float>   *elePFPhoIso_lead;
   vector<float>   *elePFNeuIso_lead;
   vector<float>   *elePFPUIso_lead;
   vector<float>   *elePFClusEcalIso_lead;
   vector<float>   *elePFClusHcalIso_lead;
   vector<float>   *eleIDMVAIso_lead;
   vector<float>   *eleIDMVANoIso_lead;
   vector<float>   *eledEtaseedAtVtx_lead;
   vector<float>   *eleE1x5_lead;
   vector<float>   *eleE2x5_lead;
   vector<float>   *eleE5x5_lead;
   vector<float>   *eleE1x5Full5x5_lead;
   vector<float>   *eleE2x5Full5x5_lead;
   vector<float>   *eleE5x5Full5x5_lead;
   vector<float>   *eleR9Full5x5_lead;
   vector<int>     *eleEcalDrivenSeed_lead;
   vector<float>   *eleDr03EcalRecHitSumEt_lead;
   vector<float>   *eleDr03HcalDepth1TowerSumEt_lead;
   vector<float>   *eleDr03HcalDepth2TowerSumEt_lead;
   vector<float>   *eleDr03HcalTowerSumEt_lead;
   vector<float>   *eleDr03TkSumPt_lead;
   vector<float>   *elecaloEnergy_lead;
   vector<float>   *eleTrkdxy_lead;
   vector<float>   *eleKFHits_lead;
   vector<float>   *eleKFChi2_lead;
   vector<float>   *eleGSFChi2_lead;
   vector<unsigned short> *eleIDbit_lead;
   vector<int>     *eleCharge_sublead;
   vector<int>     *eleChargeConsistent_sublead;
   vector<float>   *eleEn_sublead;
   vector<float>   *eleSCEn_sublead;
   vector<float>   *eleEcalEn_sublead;
   vector<float>   *eleESEnP1_sublead;
   vector<float>   *eleESEnP2_sublead;
   vector<float>   *eleD0_sublead;
   vector<float>   *eleDz_sublead;
   vector<float>   *eleD0Error_sublead;
   vector<float>   *eleDzError_sublead;
   vector<float>   *eleSIP_sublead;
   vector<float>   *elePt_sublead;
   vector<float>   *eleEta_sublead;
   vector<float>   *elePhi_sublead;
   vector<float>   *eleR9_sublead;
   vector<float>   *eleCalibPt_sublead;
   vector<float>   *eleCalibEn_sublead;
   vector<float>   *eleSCEta_sublead;
   vector<float>   *eleSCPhi_sublead;
   vector<float>   *eleSCRawEn_sublead;
   vector<float>   *eleSCEtaWidth_sublead;
   vector<float>   *eleSCPhiWidth_sublead;
   vector<float>   *eleHoverE_sublead;
   vector<float>   *eleEoverP_sublead;
   vector<float>   *eleEoverPout_sublead;
   vector<float>   *eleEoverPInv_sublead;
   vector<float>   *eleBrem_sublead;
   vector<float>   *eledEtaAtVtx_sublead;
   vector<float>   *eledPhiAtVtx_sublead;
   vector<float>   *eledEtaAtCalo_sublead;
   vector<float>   *eleSigmaIEtaIEtaFull5x5_sublead;
   vector<float>   *eleSigmaIPhiIPhiFull5x5_sublead;
   vector<int>     *eleConvVeto_sublead;
   vector<int>     *eleMissHits_sublead;
   vector<float>   *eleESEffSigmaRR_sublead;
   vector<float>   *elePFChIso_sublead;
   vector<float>   *elePFPhoIso_sublead;
   vector<float>   *elePFNeuIso_sublead;
   vector<float>   *elePFPUIso_sublead;
   vector<float>   *elePFClusEcalIso_sublead;
   vector<float>   *elePFClusHcalIso_sublead;
   vector<float>   *eleIDMVAIso_sublead;
   vector<float>   *eleIDMVANoIso_sublead;
   vector<float>   *eledEtaseedAtVtx_sublead;
   vector<float>   *eleE1x5_sublead;
   vector<float>   *eleE2x5_sublead;
   vector<float>   *eleE5x5_sublead;
   vector<float>   *eleE1x5Full5x5_sublead;
   vector<float>   *eleE2x5Full5x5_sublead;
   vector<float>   *eleE5x5Full5x5_sublead;
   vector<float>   *eleR9Full5x5_sublead;
   vector<int>     *eleEcalDrivenSeed_sublead;
   vector<float>   *eleDr03EcalRecHitSumEt_sublead;
   vector<float>   *eleDr03HcalDepth1TowerSumEt_sublead;
   vector<float>   *eleDr03HcalDepth2TowerSumEt_sublead;
   vector<float>   *eleDr03HcalTowerSumEt_sublead;
   vector<float>   *eleDr03TkSumPt_sublead;
   vector<float>   *elecaloEnergy_sublead;
   vector<float>   *eleTrkdxy_sublead;
   vector<float>   *eleKFHits_sublead;
   vector<float>   *eleKFChi2_sublead;
   vector<float>   *eleGSFChi2_sublead;
   vector<unsigned short> *eleIDbit_sublead;
   vector<float>   *eleSvChi2;
   vector<float>   *eleSvNDOF;
   vector<float>   *eleSvProb;
   vector<float>   *eleSvX;
   vector<float>   *eleSvY;
   vector<float>   *eleSvZ;
   vector<float>   *eleSvXError;
   vector<float>   *eleSvYError;
   vector<float>   *eleSvZError;
   vector<float>   *eleSvMass;
   vector<float>   *eleSvCtxy;
   vector<float>   *eleSvCosAngle;
   vector<float>   *eleSvLxy;
   vector<float>   *eleSvLxyError;
   vector<int>     *kaonEECharge_lead;
   vector<float>   *kaonEED0_lead;
   vector<float>   *kaonEEDz_lead;
   vector<float>   *kaonEED0Error_lead;
   vector<float>   *kaonEEDzError_lead;
   vector<float>   *kaonEEPt_lead;
   vector<float>   *kaonEEEta_lead;
   vector<float>   *kaonEEPhi_lead;
   vector<float>   *kaonEEVx_lead;
   vector<float>   *kaonEEVy_lead;
   vector<float>   *kaonEEVz_lead;
   vector<float>   *kaonEETrkChi2_lead;
   vector<float>   *kaonEETrkNDOF_lead;
   vector<float>   *kaonEETrkNormChi2_lead;
   vector<int>     *kaonEECharge_sublead;
   vector<float>   *kaonEED0_sublead;
   vector<float>   *kaonEEDz_sublead;
   vector<float>   *kaonEED0Error_sublead;
   vector<float>   *kaonEEDzError_sublead;
   vector<float>   *kaonEEPt_sublead;
   vector<float>   *kaonEEEta_sublead;
   vector<float>   *kaonEEPhi_sublead;
   vector<float>   *kaonEEVx_sublead;
   vector<float>   *kaonEEVy_sublead;
   vector<float>   *kaonEEVz_sublead;
   vector<float>   *kaonEETrkChi2_sublead;
   vector<float>   *kaonEETrkNDOF_sublead;
   vector<float>   *kaonEETrkNormChi2_sublead;
   vector<float>   *bsEEdRele;
   vector<float>   *bsEEdRkaon;
   vector<float>   *bsEEdRJpsiPhi;
   vector<float>   *bsEEJpsiMass;
   vector<float>   *bsEEPhiMass;
   vector<float>   *bsEEBsMass;
   Int_t           nMu;
   Bool_t          matchedTrg;
   vector<float>   *muPt_lead;
   vector<float>   *muEn_lead;
   vector<float>   *muEta_lead;
   vector<float>   *muPhi_lead;
   vector<int>     *muCharge_lead;
   vector<int>     *muType_lead;
   vector<unsigned short> *muIDbit_lead;
   vector<float>   *muD0_lead;
   vector<float>   *muDz_lead;
   vector<float>   *muSIP_lead;
   vector<float>   *muD0Error_lead;
   vector<float>   *muDzError_lead;
   vector<float>   *muChi2NDF_lead;
   vector<float>   *muInnerD0_lead;
   vector<float>   *muInnerDz_lead;
   vector<int>     *muTrkLayers_lead;
   vector<int>     *muPixelLayers_lead;
   vector<int>     *muPixelHits_lead;
   vector<int>     *muMuonHits_lead;
   vector<int>     *muStations_lead;
   vector<int>     *muMatches_lead;
   vector<int>     *muTrkQuality_lead;
   vector<float>   *muIsoTrk_lead;
   vector<float>   *muPFChIso_lead;
   vector<float>   *muPFPhoIso_lead;
   vector<float>   *muPFNeuIso_lead;
   vector<float>   *muPFPUIso_lead;
   vector<float>   *muPFChIso03_lead;
   vector<float>   *muPFPhoIso03_lead;
   vector<float>   *muPFNeuIso03_lead;
   vector<float>   *muPFPUIso03_lead;
   vector<bool>    *muFiredTrgs_lead;
   vector<ULong64_t> *muFiredL1Trgs_lead;
   vector<float>   *muInnervalidFraction_lead;
   vector<float>   *musegmentCompatibility_lead;
   vector<float>   *muchi2LocalPosition_lead;
   vector<float>   *mutrkKink_lead;
   vector<float>   *muBestTrkPtError_lead;
   vector<float>   *muBestTrkPt_lead;
   vector<int>     *muBestTrkType_lead;
   vector<float>   *muTrkNormChi2_lead;
   vector<float>   *muIDPatMVA_lead;
   vector<float>   *muIDPatSoftMVA_lead;
   vector<unsigned short> *muIDSelectorBit_lead;
   vector<float>   *muPt_sublead;
   vector<float>   *muEn_sublead;
   vector<float>   *muEta_sublead;
   vector<float>   *muPhi_sublead;
   vector<int>     *muCharge_sublead;
   vector<int>     *muType_sublead;
   vector<unsigned short> *muIDbit_sublead;
   vector<float>   *muD0_sublead;
   vector<float>   *muDz_sublead;
   vector<float>   *muSIP_sublead;
   vector<float>   *muD0Error_sublead;
   vector<float>   *muDzError_sublead;
   vector<float>   *muChi2NDF_sublead;
   vector<float>   *muInnerD0_sublead;
   vector<float>   *muInnerDz_sublead;
   vector<int>     *muTrkLayers_sublead;
   vector<int>     *muPixelLayers_sublead;
   vector<int>     *muPixelHits_sublead;
   vector<int>     *muMuonHits_sublead;
   vector<int>     *muStations_sublead;
   vector<int>     *muMatches_sublead;
   vector<int>     *muTrkQuality_sublead;
   vector<float>   *muIsoTrk_sublead;
   vector<float>   *muPFChIso_sublead;
   vector<float>   *muPFPhoIso_sublead;
   vector<float>   *muPFNeuIso_sublead;
   vector<float>   *muPFPUIso_sublead;
   vector<float>   *muPFChIso03_sublead;
   vector<float>   *muPFPhoIso03_sublead;
   vector<float>   *muPFNeuIso03_sublead;
   vector<float>   *muPFPUIso03_sublead;
   vector<bool>    *muFiredTrgs_sublead;
   vector<ULong64_t> *muFiredL1Trgs_sublead;
   vector<float>   *muInnervalidFraction_sublead;
   vector<float>   *musegmentCompatibility_sublead;
   vector<float>   *muchi2LocalPosition_sublead;
   vector<float>   *mutrkKink_sublead;
   vector<float>   *muBestTrkPtError_sublead;
   vector<float>   *muBestTrkPt_sublead;
   vector<int>     *muBestTrkType_sublead;
   vector<float>   *muTrkNormChi2_sublead;
   vector<float>   *muIDPatMVA_sublead;
   vector<float>   *muIDPatSoftMVA_sublead;
   vector<unsigned short> *muIDSelectorBit_sublead;
   vector<float>   *muSvChi2;
   vector<float>   *muSvNDOF;
   vector<float>   *muSvX;
   vector<float>   *muSvY;
   vector<float>   *muSvZ;
   vector<float>   *muSvXError;
   vector<float>   *muSvYError;
   vector<float>   *muSvZError;
   vector<float>   *muSvMass;
   vector<float>   *muSvCtxy;
   vector<float>   *muSvCosAngle;
   vector<int>     *kaonMMCharge_lead;
   vector<float>   *kaonMMD0_lead;
   vector<float>   *kaonMMDz_lead;
   vector<float>   *kaonMMD0Error_lead;
   vector<float>   *kaonMMDzError_lead;
   vector<float>   *kaonMMPt_lead;
   vector<float>   *kaonMMEta_lead;
   vector<float>   *kaonMMPhi_lead;
   vector<float>   *kaonMMVx_lead;
   vector<float>   *kaonMMVy_lead;
   vector<float>   *kaonMMVz_lead;
   vector<float>   *kaonMMEn_lead;
   vector<float>   *kaonMMTrkChi2_lead;
   vector<float>   *kaonMMTrkNDOF_lead;
   vector<float>   *kaonMMTrkNormChi2_lead;
   vector<int>     *kaonMMCharge_sublead;
   vector<float>   *kaonMMD0_sublead;
   vector<float>   *kaonMMDz_sublead;
   vector<float>   *kaonMMD0Error_sublead;
   vector<float>   *kaonMMDzError_sublead;
   vector<float>   *kaonMMPt_sublead;
   vector<float>   *kaonMMEta_sublead;
   vector<float>   *kaonMMPhi_sublead;
   vector<float>   *kaonMMVx_sublead;
   vector<float>   *kaonMMVy_sublead;
   vector<float>   *kaonMMVz_sublead;
   vector<float>   *kaonMMEn_sublead;
   vector<float>   *kaonMMTrkChi2_sublead;
   vector<float>   *kaonMMTrkNDOF_sublead;
   vector<float>   *kaonMMTrkNormChi2_sublead;
   vector<float>   *bsMMdRmu;
   vector<float>   *bsMMdRkaon;
   vector<float>   *bsMMdRJpsiPhi;
   vector<float>   *bsMMJpsiMass;
   vector<float>   *bsMMPhiMass;
   vector<float>   *bsMMBsMass;

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
   //TBranch        *b_npfPho;   //!
   //TBranch        *b_pfphoEt;   //!
   //TBranch        *b_pfphoEta;   //!
   //TBranch        *b_pfphoPhi;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge_lead;   //!
   TBranch        *b_eleChargeConsistent_lead;   //!
   TBranch        *b_eleEn_lead;   //!
   TBranch        *b_eleSCEn_lead;   //!
   TBranch        *b_eleEcalEn_lead;   //!
   TBranch        *b_eleESEnP1_lead;   //!
   TBranch        *b_eleESEnP2_lead;   //!
   TBranch        *b_eleD0_lead;   //!
   TBranch        *b_eleDz_lead;   //!
   TBranch        *b_eleD0Error_lead;   //!
   TBranch        *b_eleDzError_lead;   //!
   TBranch        *b_eleSIP_lead;   //!
   TBranch        *b_elePt_lead;   //!
   TBranch        *b_eleEta_lead;   //!
   TBranch        *b_elePhi_lead;   //!
   TBranch        *b_eleR9_lead;   //!
   TBranch        *b_eleCalibPt_lead;   //!
   TBranch        *b_eleCalibEn_lead;   //!
   TBranch        *b_eleSCEta_lead;   //!
   TBranch        *b_eleSCPhi_lead;   //!
   TBranch        *b_eleSCRawEn_lead;   //!
   TBranch        *b_eleSCEtaWidth_lead;   //!
   TBranch        *b_eleSCPhiWidth_lead;   //!
   TBranch        *b_eleHoverE_lead;   //!
   TBranch        *b_eleEoverP_lead;   //!
   TBranch        *b_eleEoverPout_lead;   //!
   TBranch        *b_eleEoverPInv_lead;   //!
   TBranch        *b_eleBrem_lead;   //!
   TBranch        *b_eledEtaAtVtx_lead;   //!
   TBranch        *b_eledPhiAtVtx_lead;   //!
   TBranch        *b_eledEtaAtCalo_lead;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5_lead;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5_lead;   //!
   TBranch        *b_eleConvVeto_lead;   //!
   TBranch        *b_eleMissHits_lead;   //!
   TBranch        *b_eleESEffSigmaRR_lead;   //!
   TBranch        *b_elePFChIso_lead;   //!
   TBranch        *b_elePFPhoIso_lead;   //!
   TBranch        *b_elePFNeuIso_lead;   //!
   TBranch        *b_elePFPUIso_lead;   //!
   TBranch        *b_elePFClusEcalIso_lead;   //!
   TBranch        *b_elePFClusHcalIso_lead;   //!
   TBranch        *b_eleIDMVAIso_lead;   //!
   TBranch        *b_eleIDMVANoIso_lead;   //!
   TBranch        *b_eledEtaseedAtVtx_lead;   //!
   TBranch        *b_eleE1x5_lead;   //!
   TBranch        *b_eleE2x5_lead;   //!
   TBranch        *b_eleE5x5_lead;   //!
   TBranch        *b_eleE1x5Full5x5_lead;   //!
   TBranch        *b_eleE2x5Full5x5_lead;   //!
   TBranch        *b_eleE5x5Full5x5_lead;   //!
   TBranch        *b_eleR9Full5x5_lead;   //!
   TBranch        *b_eleEcalDrivenSeed_lead;   //!
   TBranch        *b_eleDr03EcalRecHitSumEt_lead;   //!
   TBranch        *b_eleDr03HcalDepth1TowerSumEt_lead;   //!
   TBranch        *b_eleDr03HcalDepth2TowerSumEt_lead;   //!
   TBranch        *b_eleDr03HcalTowerSumEt_lead;   //!
   TBranch        *b_eleDr03TkSumPt_lead;   //!
   TBranch        *b_elecaloEnergy_lead;   //!
   TBranch        *b_eleTrkdxy_lead;   //!
   TBranch        *b_eleKFHits_lead;   //!
   TBranch        *b_eleKFChi2_lead;   //!
   TBranch        *b_eleGSFChi2_lead;   //!
   TBranch        *b_eleIDbit_lead;   //!
   TBranch        *b_eleCharge_sublead;   //!
   TBranch        *b_eleChargeConsistent_sublead;   //!
   TBranch        *b_eleEn_sublead;   //!
   TBranch        *b_eleSCEn_sublead;   //!
   TBranch        *b_eleEcalEn_sublead;   //!
   TBranch        *b_eleESEnP1_sublead;   //!
   TBranch        *b_eleESEnP2_sublead;   //!
   TBranch        *b_eleD0_sublead;   //!
   TBranch        *b_eleDz_sublead;   //!
   TBranch        *b_eleD0Error_sublead;   //!
   TBranch        *b_eleDzError_sublead;   //!
   TBranch        *b_eleSIP_sublead;   //!
   TBranch        *b_elePt_sublead;   //!
   TBranch        *b_eleEta_sublead;   //!
   TBranch        *b_elePhi_sublead;   //!
   TBranch        *b_eleR9_sublead;   //!
   TBranch        *b_eleCalibPt_sublead;   //!
   TBranch        *b_eleCalibEn_sublead;   //!
   TBranch        *b_eleSCEta_sublead;   //!
   TBranch        *b_eleSCPhi_sublead;   //!
   TBranch        *b_eleSCRawEn_sublead;   //!
   TBranch        *b_eleSCEtaWidth_sublead;   //!
   TBranch        *b_eleSCPhiWidth_sublead;   //!
   TBranch        *b_eleHoverE_sublead;   //!
   TBranch        *b_eleEoverP_sublead;   //!
   TBranch        *b_eleEoverPout_sublead;   //!
   TBranch        *b_eleEoverPInv_sublead;   //!
   TBranch        *b_eleBrem_sublead;   //!
   TBranch        *b_eledEtaAtVtx_sublead;   //!
   TBranch        *b_eledPhiAtVtx_sublead;   //!
   TBranch        *b_eledEtaAtCalo_sublead;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5_sublead;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5_sublead;   //!
   TBranch        *b_eleConvVeto_sublead;   //!
   TBranch        *b_eleMissHits_sublead;   //!
   TBranch        *b_eleESEffSigmaRR_sublead;   //!
   TBranch        *b_elePFChIso_sublead;   //!
   TBranch        *b_elePFPhoIso_sublead;   //!
   TBranch        *b_elePFNeuIso_sublead;   //!
   TBranch        *b_elePFPUIso_sublead;   //!
   TBranch        *b_elePFClusEcalIso_sublead;   //!
   TBranch        *b_elePFClusHcalIso_sublead;   //!
   TBranch        *b_eleIDMVAIso_sublead;   //!
   TBranch        *b_eleIDMVANoIso_sublead;   //!
   TBranch        *b_eledEtaseedAtVtx_sublead;   //!
   TBranch        *b_eleE1x5_sublead;   //!
   TBranch        *b_eleE2x5_sublead;   //!
   TBranch        *b_eleE5x5_sublead;   //!
   TBranch        *b_eleE1x5Full5x5_sublead;   //!
   TBranch        *b_eleE2x5Full5x5_sublead;   //!
   TBranch        *b_eleE5x5Full5x5_sublead;   //!
   TBranch        *b_eleR9Full5x5_sublead;   //!
   TBranch        *b_eleEcalDrivenSeed_sublead;   //!
   TBranch        *b_eleDr03EcalRecHitSumEt_sublead;   //!
   TBranch        *b_eleDr03HcalDepth1TowerSumEt_sublead;   //!
   TBranch        *b_eleDr03HcalDepth2TowerSumEt_sublead;   //!
   TBranch        *b_eleDr03HcalTowerSumEt_sublead;   //!
   TBranch        *b_eleDr03TkSumPt_sublead;   //!
   TBranch        *b_elecaloEnergy_sublead;   //!
   TBranch        *b_eleTrkdxy_sublead;   //!
   TBranch        *b_eleKFHits_sublead;   //!
   TBranch        *b_eleKFChi2_sublead;   //!
   TBranch        *b_eleGSFChi2_sublead;   //!
   TBranch        *b_eleIDbit_sublead;   //!
   TBranch        *b_eleSvChi2;   //!
   TBranch        *b_eleSvNDOF;   //!
   TBranch        *b_eleSvProb;   //!
   TBranch        *b_eleSvX;   //!
   TBranch        *b_eleSvY;   //!
   TBranch        *b_eleSvZ;   //!
   TBranch        *b_eleSvXError;   //!
   TBranch        *b_eleSvYError;   //!
   TBranch        *b_eleSvZError;   //!
   TBranch        *b_eleSvMass;   //!
   TBranch        *b_eleSvCtxy;   //!
   TBranch        *b_eleSvCosAngle;   //!
   TBranch        *b_eleSvLxy;   //!
   TBranch        *b_eleSvLxyError;   //!
   TBranch        *b_kaonEECharge_lead;   //!
   TBranch        *b_kaonEED0_lead;   //!
   TBranch        *b_kaonEEDz_lead;   //!
   TBranch        *b_kaonEED0Error_lead;   //!
   TBranch        *b_kaonEEDzError_lead;   //!
   TBranch        *b_kaonEEPt_lead;   //!
   TBranch        *b_kaonEEEta_lead;   //!
   TBranch        *b_kaonEEPhi_lead;   //!
   TBranch        *b_kaonEEVx_lead;   //!
   TBranch        *b_kaonEEVy_lead;   //!
   TBranch        *b_kaonEEVz_lead;   //!
   TBranch        *b_kaonEETrkChi2_lead;   //!
   TBranch        *b_kaonEETrkNDOF_lead;   //!
   TBranch        *b_kaonEETrkNormChi2_lead;   //!
   TBranch        *b_kaonEECharge_sublead;   //!
   TBranch        *b_kaonEED0_sublead;   //!
   TBranch        *b_kaonEEDz_sublead;   //!
   TBranch        *b_kaonEED0Error_sublead;   //!
   TBranch        *b_kaonEEDzError_sublead;   //!
   TBranch        *b_kaonEEPt_sublead;   //!
   TBranch        *b_kaonEEEta_sublead;   //!
   TBranch        *b_kaonEEPhi_sublead;   //!
   TBranch        *b_kaonEEVx_sublead;   //!
   TBranch        *b_kaonEEVy_sublead;   //!
   TBranch        *b_kaonEEVz_sublead;   //!
   TBranch        *b_kaonEETrkChi2_sublead;   //!
   TBranch        *b_kaonEETrkNDOF_sublead;   //!
   TBranch        *b_kaonEETrkNormChi2_sublead;   //!
   TBranch        *b_bsEEdRele;   //!
   TBranch        *b_bsEEdRkaon;   //!
   TBranch        *b_bsEEdRJpsiPhi;   //!
   TBranch        *b_bsEEJpsiMass;   //!
   TBranch        *b_bsEEPhiMass;   //!
   TBranch        *b_bsEEBsMass;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_matchedTrg;   //!
   TBranch        *b_muPt_lead;   //!
   TBranch        *b_muEn_lead;   //!
   TBranch        *b_muEta_lead;   //!
   TBranch        *b_muPhi_lead;   //!
   TBranch        *b_muCharge_lead;   //!
   TBranch        *b_muType_lead;   //!
   TBranch        *b_muIDbit_lead;   //!
   TBranch        *b_muD0_lead;   //!
   TBranch        *b_muDz_lead;   //!
   TBranch        *b_muSIP_lead;   //!
   TBranch        *b_muD0Error_lead;   //!
   TBranch        *b_muDzError_lead;   //!
   TBranch        *b_muChi2NDF_lead;   //!
   TBranch        *b_muInnerD0_lead;   //!
   TBranch        *b_muInnerDz_lead;   //!
   TBranch        *b_muTrkLayers_lead;   //!
   TBranch        *b_muPixelLayers_lead;   //!
   TBranch        *b_muPixelHits_lead;   //!
   TBranch        *b_muMuonHits_lead;   //!
   TBranch        *b_muStations_lead;   //!
   TBranch        *b_muMatches_lead;   //!
   TBranch        *b_muTrkQuality_lead;   //!
   TBranch        *b_muIsoTrk_lead;   //!
   TBranch        *b_muPFChIso_lead;   //!
   TBranch        *b_muPFPhoIso_lead;   //!
   TBranch        *b_muPFNeuIso_lead;   //!
   TBranch        *b_muPFPUIso_lead;   //!
   TBranch        *b_muPFChIso03_lead;   //!
   TBranch        *b_muPFPhoIso03_lead;   //!
   TBranch        *b_muPFNeuIso03_lead;   //!
   TBranch        *b_muPFPUIso03_lead;   //!
   TBranch        *b_muFiredTrgs_lead;   //!
   TBranch        *b_muFiredL1Trgs_lead;   //!
   TBranch        *b_muInnervalidFraction_lead;   //!
   TBranch        *b_musegmentCompatibility_lead;   //!
   TBranch        *b_muchi2LocalPosition_lead;   //!
   TBranch        *b_mutrkKink_lead;   //!
   TBranch        *b_muBestTrkPtError_lead;   //!
   TBranch        *b_muBestTrkPt_lead;   //!
   TBranch        *b_muBestTrkType_lead;   //!
   TBranch        *b_muTrkNormChi2_lead;   //!
   TBranch        *b_muIDPatMVA_lead;   //!
   TBranch        *b_muIDPatSoftMVA_lead;   //!
   TBranch        *b_muIDSelectorBit_lead;   //!
   TBranch        *b_muPt_sublead;   //!
   TBranch        *b_muEn_sublead;   //!
   TBranch        *b_muEta_sublead;   //!
   TBranch        *b_muPhi_sublead;   //!
   TBranch        *b_muCharge_sublead;   //!
   TBranch        *b_muType_sublead;   //!
   TBranch        *b_muIDbit_sublead;   //!
   TBranch        *b_muD0_sublead;   //!
   TBranch        *b_muDz_sublead;   //!
   TBranch        *b_muSIP_sublead;   //!
   TBranch        *b_muD0Error_sublead;   //!
   TBranch        *b_muDzError_sublead;   //!
   TBranch        *b_muChi2NDF_sublead;   //!
   TBranch        *b_muInnerD0_sublead;   //!
   TBranch        *b_muInnerDz_sublead;   //!
   TBranch        *b_muTrkLayers_sublead;   //!
   TBranch        *b_muPixelLayers_sublead;   //!
   TBranch        *b_muPixelHits_sublead;   //!
   TBranch        *b_muMuonHits_sublead;   //!
   TBranch        *b_muStations_sublead;   //!
   TBranch        *b_muMatches_sublead;   //!
   TBranch        *b_muTrkQuality_sublead;   //!
   TBranch        *b_muIsoTrk_sublead;   //!
   TBranch        *b_muPFChIso_sublead;   //!
   TBranch        *b_muPFPhoIso_sublead;   //!
   TBranch        *b_muPFNeuIso_sublead;   //!
   TBranch        *b_muPFPUIso_sublead;   //!
   TBranch        *b_muPFChIso03_sublead;   //!
   TBranch        *b_muPFPhoIso03_sublead;   //!
   TBranch        *b_muPFNeuIso03_sublead;   //!
   TBranch        *b_muPFPUIso03_sublead;   //!
   TBranch        *b_muFiredTrgs_sublead;   //!
   TBranch        *b_muFiredL1Trgs_sublead;   //!
   TBranch        *b_muInnervalidFraction_sublead;   //!
   TBranch        *b_musegmentCompatibility_sublead;   //!
   TBranch        *b_muchi2LocalPosition_sublead;   //!
   TBranch        *b_mutrkKink_sublead;   //!
   TBranch        *b_muBestTrkPtError_sublead;   //!
   TBranch        *b_muBestTrkPt_sublead;   //!
   TBranch        *b_muBestTrkType_sublead;   //!
   TBranch        *b_muTrkNormChi2_sublead;   //!
   TBranch        *b_muIDPatMVA_sublead;   //!
   TBranch        *b_muIDPatSoftMVA_sublead;   //!
   TBranch        *b_muIDSelectorBit_sublead;   //!
   TBranch        *b_muSvChi2;   //!
   TBranch        *b_muSvNDOF;   //!
   TBranch        *b_muSvX;   //!
   TBranch        *b_muSvY;   //!
   TBranch        *b_muSvZ;   //!
   TBranch        *b_muSvXError;   //!
   TBranch        *b_muSvYError;   //!
   TBranch        *b_muSvZError;   //!
   TBranch        *b_muSvMass;   //!
   TBranch        *b_muSvCtxy;   //!
   TBranch        *b_muSvCosAngle;   //!
   TBranch        *b_kaonMMCharge_lead;   //!
   TBranch        *b_kaonMMD0_lead;   //!
   TBranch        *b_kaonMMDz_lead;   //!
   TBranch        *b_kaonMMD0Error_lead;   //!
   TBranch        *b_kaonMMDzError_lead;   //!
   TBranch        *b_kaonMMPt_lead;   //!
   TBranch        *b_kaonMMEta_lead;   //!
   TBranch        *b_kaonMMPhi_lead;   //!
   TBranch        *b_kaonMMVx_lead;   //!
   TBranch        *b_kaonMMVy_lead;   //!
   TBranch        *b_kaonMMVz_lead;   //!
   TBranch        *b_kaonMMEn_lead;   //!
   TBranch        *b_kaonMMTrkChi2_lead;   //!
   TBranch        *b_kaonMMTrkNDOF_lead;   //!
   TBranch        *b_kaonMMTrkNormChi2_lead;   //!
   TBranch        *b_kaonMMCharge_sublead;   //!
   TBranch        *b_kaonMMD0_sublead;   //!
   TBranch        *b_kaonMMDz_sublead;   //!
   TBranch        *b_kaonMMD0Error_sublead;   //!
   TBranch        *b_kaonMMDzError_sublead;   //!
   TBranch        *b_kaonMMPt_sublead;   //!
   TBranch        *b_kaonMMEta_sublead;   //!
   TBranch        *b_kaonMMPhi_sublead;   //!
   TBranch        *b_kaonMMVx_sublead;   //!
   TBranch        *b_kaonMMVy_sublead;   //!
   TBranch        *b_kaonMMVz_sublead;   //!
   TBranch        *b_kaonMMEn_sublead;   //!
   TBranch        *b_kaonMMTrkChi2_sublead;   //!
   TBranch        *b_kaonMMTrkNDOF_sublead;   //!
   TBranch        *b_kaonMMTrkNormChi2_sublead;   //!
   TBranch        *b_bsMMdRmu;   //!
   TBranch        *b_bsMMdRkaon;   //!
   TBranch        *b_bsMMdRJpsiPhi;   //!
   TBranch        *b_bsMMJpsiMass;   //!
   TBranch        *b_bsMMPhiMass;   //!
   TBranch        *b_bsMMBsMass;   //!

   Float_t elecM = 0.000510998;
   Float_t muonM = 0.1056583745;
   Float_t kaonM = 0.493677;
   Float_t bsM = 5.3663;
   Float_t jpsilow = 2.9;
   Float_t jpsiup = 3.45;
   Float_t philow = 1.01;
   Float_t phiup = 1.04;
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

   BsPhiLLTupleTree(TTree *tree=0);
   virtual ~BsPhiLLTupleTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outputfile, Int_t maxevents, Float_t mvaCut);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Float_t          deltaPhi(Float_t phi1, Float_t phi2);
};

#endif

#ifdef BsPhiLLTupleTree_cxx
BsPhiLLTupleTree::BsPhiLLTupleTree(TTree *tree) : fChain(0) 
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

BsPhiLLTupleTree::~BsPhiLLTupleTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BsPhiLLTupleTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BsPhiLLTupleTree::LoadTree(Long64_t entry)
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

void BsPhiLLTupleTree::Init(TTree *tree)
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
   //pfphoEt = 0;
   //pfphoEta = 0;
   //pfphoPhi = 0;
   eleCharge_lead = 0;
   eleChargeConsistent_lead = 0;
   eleEn_lead = 0;
   eleSCEn_lead = 0;
   eleEcalEn_lead = 0;
   eleESEnP1_lead = 0;
   eleESEnP2_lead = 0;
   eleD0_lead = 0;
   eleDz_lead = 0;
   eleD0Error_lead = 0;
   eleDzError_lead = 0;
   eleSIP_lead = 0;
   elePt_lead = 0;
   eleEta_lead = 0;
   elePhi_lead = 0;
   eleR9_lead = 0;
   eleCalibPt_lead = 0;
   eleCalibEn_lead = 0;
   eleSCEta_lead = 0;
   eleSCPhi_lead = 0;
   eleSCRawEn_lead = 0;
   eleSCEtaWidth_lead = 0;
   eleSCPhiWidth_lead = 0;
   eleHoverE_lead = 0;
   eleEoverP_lead = 0;
   eleEoverPout_lead = 0;
   eleEoverPInv_lead = 0;
   eleBrem_lead = 0;
   eledEtaAtVtx_lead = 0;
   eledPhiAtVtx_lead = 0;
   eledEtaAtCalo_lead = 0;
   eleSigmaIEtaIEtaFull5x5_lead = 0;
   eleSigmaIPhiIPhiFull5x5_lead = 0;
   eleConvVeto_lead = 0;
   eleMissHits_lead = 0;
   eleESEffSigmaRR_lead = 0;
   elePFChIso_lead = 0;
   elePFPhoIso_lead = 0;
   elePFNeuIso_lead = 0;
   elePFPUIso_lead = 0;
   elePFClusEcalIso_lead = 0;
   elePFClusHcalIso_lead = 0;
   eleIDMVAIso_lead = 0;
   eleIDMVANoIso_lead = 0;
   eledEtaseedAtVtx_lead = 0;
   eleE1x5_lead = 0;
   eleE2x5_lead = 0;
   eleE5x5_lead = 0;
   eleE1x5Full5x5_lead = 0;
   eleE2x5Full5x5_lead = 0;
   eleE5x5Full5x5_lead = 0;
   eleR9Full5x5_lead = 0;
   eleEcalDrivenSeed_lead = 0;
   eleDr03EcalRecHitSumEt_lead = 0;
   eleDr03HcalDepth1TowerSumEt_lead = 0;
   eleDr03HcalDepth2TowerSumEt_lead = 0;
   eleDr03HcalTowerSumEt_lead = 0;
   eleDr03TkSumPt_lead = 0;
   elecaloEnergy_lead = 0;
   eleTrkdxy_lead = 0;
   eleKFHits_lead = 0;
   eleKFChi2_lead = 0;
   eleGSFChi2_lead = 0;
   eleIDbit_lead = 0;
   eleCharge_sublead = 0;
   eleChargeConsistent_sublead = 0;
   eleEn_sublead = 0;
   eleSCEn_sublead = 0;
   eleEcalEn_sublead = 0;
   eleESEnP1_sublead = 0;
   eleESEnP2_sublead = 0;
   eleD0_sublead = 0;
   eleDz_sublead = 0;
   eleD0Error_sublead = 0;
   eleDzError_sublead = 0;
   eleSIP_sublead = 0;
   elePt_sublead = 0;
   eleEta_sublead = 0;
   elePhi_sublead = 0;
   eleR9_sublead = 0;
   eleCalibPt_sublead = 0;
   eleCalibEn_sublead = 0;
   eleSCEta_sublead = 0;
   eleSCPhi_sublead = 0;
   eleSCRawEn_sublead = 0;
   eleSCEtaWidth_sublead = 0;
   eleSCPhiWidth_sublead = 0;
   eleHoverE_sublead = 0;
   eleEoverP_sublead = 0;
   eleEoverPout_sublead = 0;
   eleEoverPInv_sublead = 0;
   eleBrem_sublead = 0;
   eledEtaAtVtx_sublead = 0;
   eledPhiAtVtx_sublead = 0;
   eledEtaAtCalo_sublead = 0;
   eleSigmaIEtaIEtaFull5x5_sublead = 0;
   eleSigmaIPhiIPhiFull5x5_sublead = 0;
   eleConvVeto_sublead = 0;
   eleMissHits_sublead = 0;
   eleESEffSigmaRR_sublead = 0;
   elePFChIso_sublead = 0;
   elePFPhoIso_sublead = 0;
   elePFNeuIso_sublead = 0;
   elePFPUIso_sublead = 0;
   elePFClusEcalIso_sublead = 0;
   elePFClusHcalIso_sublead = 0;
   eleIDMVAIso_sublead = 0;
   eleIDMVANoIso_sublead = 0;
   eledEtaseedAtVtx_sublead = 0;
   eleE1x5_sublead = 0;
   eleE2x5_sublead = 0;
   eleE5x5_sublead = 0;
   eleE1x5Full5x5_sublead = 0;
   eleE2x5Full5x5_sublead = 0;
   eleE5x5Full5x5_sublead = 0;
   eleR9Full5x5_sublead = 0;
   eleEcalDrivenSeed_sublead = 0;
   eleDr03EcalRecHitSumEt_sublead = 0;
   eleDr03HcalDepth1TowerSumEt_sublead = 0;
   eleDr03HcalDepth2TowerSumEt_sublead = 0;
   eleDr03HcalTowerSumEt_sublead = 0;
   eleDr03TkSumPt_sublead = 0;
   elecaloEnergy_sublead = 0;
   eleTrkdxy_sublead = 0;
   eleKFHits_sublead = 0;
   eleKFChi2_sublead = 0;
   eleGSFChi2_sublead = 0;
   eleIDbit_sublead = 0;
   eleSvChi2 = 0;
   eleSvNDOF = 0;
   eleSvProb = 0;
   eleSvX = 0;
   eleSvY = 0;
   eleSvZ = 0;
   eleSvXError = 0;
   eleSvYError = 0;
   eleSvZError = 0;
   eleSvMass = 0;
   eleSvCtxy = 0;
   eleSvCosAngle = 0;
   eleSvLxy = 0;
   eleSvLxyError = 0;
   kaonEECharge_lead = 0;
   kaonEED0_lead = 0;
   kaonEEDz_lead = 0;
   kaonEED0Error_lead = 0;
   kaonEEDzError_lead = 0;
   kaonEEPt_lead = 0;
   kaonEEEta_lead = 0;
   kaonEEPhi_lead = 0;
   kaonEEVx_lead = 0;
   kaonEEVy_lead = 0;
   kaonEEVz_lead = 0;
   kaonEETrkChi2_lead = 0;
   kaonEETrkNDOF_lead = 0;
   kaonEETrkNormChi2_lead = 0;
   kaonEECharge_sublead = 0;
   kaonEED0_sublead = 0;
   kaonEEDz_sublead = 0;
   kaonEED0Error_sublead = 0;
   kaonEEDzError_sublead = 0;
   kaonEEPt_sublead = 0;
   kaonEEEta_sublead = 0;
   kaonEEPhi_sublead = 0;
   kaonEEVx_sublead = 0;
   kaonEEVy_sublead = 0;
   kaonEEVz_sublead = 0;
   kaonEETrkChi2_sublead = 0;
   kaonEETrkNDOF_sublead = 0;
   kaonEETrkNormChi2_sublead = 0;
   bsEEdRele = 0;
   bsEEdRkaon = 0;
   bsEEdRJpsiPhi = 0;
   bsEEJpsiMass = 0;
   bsEEPhiMass = 0;
   bsEEBsMass = 0;
   muPt_lead = 0;
   muEn_lead = 0;
   muEta_lead = 0;
   muPhi_lead = 0;
   muCharge_lead = 0;
   muType_lead = 0;
   muIDbit_lead = 0;
   muD0_lead = 0;
   muDz_lead = 0;
   muSIP_lead = 0;
   muD0Error_lead = 0;
   muDzError_lead = 0;
   muChi2NDF_lead = 0;
   muInnerD0_lead = 0;
   muInnerDz_lead = 0;
   muTrkLayers_lead = 0;
   muPixelLayers_lead = 0;
   muPixelHits_lead = 0;
   muMuonHits_lead = 0;
   muStations_lead = 0;
   muMatches_lead = 0;
   muTrkQuality_lead = 0;
   muIsoTrk_lead = 0;
   muPFChIso_lead = 0;
   muPFPhoIso_lead = 0;
   muPFNeuIso_lead = 0;
   muPFPUIso_lead = 0;
   muPFChIso03_lead = 0;
   muPFPhoIso03_lead = 0;
   muPFNeuIso03_lead = 0;
   muPFPUIso03_lead = 0;
   muFiredTrgs_lead = 0;
   muFiredL1Trgs_lead = 0;
   muInnervalidFraction_lead = 0;
   musegmentCompatibility_lead = 0;
   muchi2LocalPosition_lead = 0;
   mutrkKink_lead = 0;
   muBestTrkPtError_lead = 0;
   muBestTrkPt_lead = 0;
   muBestTrkType_lead = 0;
   muTrkNormChi2_lead = 0;
   muIDPatMVA_lead = 0;
   muIDPatSoftMVA_lead = 0;
   muIDSelectorBit_lead = 0;
   muPt_sublead = 0;
   muEn_sublead = 0;
   muEta_sublead = 0;
   muPhi_sublead = 0;
   muCharge_sublead = 0;
   muType_sublead = 0;
   muIDbit_sublead = 0;
   muD0_sublead = 0;
   muDz_sublead = 0;
   muSIP_sublead = 0;
   muD0Error_sublead = 0;
   muDzError_sublead = 0;
   muChi2NDF_sublead = 0;
   muInnerD0_sublead = 0;
   muInnerDz_sublead = 0;
   muTrkLayers_sublead = 0;
   muPixelLayers_sublead = 0;
   muPixelHits_sublead = 0;
   muMuonHits_sublead = 0;
   muStations_sublead = 0;
   muMatches_sublead = 0;
   muTrkQuality_sublead = 0;
   muIsoTrk_sublead = 0;
   muPFChIso_sublead = 0;
   muPFPhoIso_sublead = 0;
   muPFNeuIso_sublead = 0;
   muPFPUIso_sublead = 0;
   muPFChIso03_sublead = 0;
   muPFPhoIso03_sublead = 0;
   muPFNeuIso03_sublead = 0;
   muPFPUIso03_sublead = 0;
   muFiredTrgs_sublead = 0;
   muFiredL1Trgs_sublead = 0;
   muInnervalidFraction_sublead = 0;
   musegmentCompatibility_sublead = 0;
   muchi2LocalPosition_sublead = 0;
   mutrkKink_sublead = 0;
   muBestTrkPtError_sublead = 0;
   muBestTrkPt_sublead = 0;
   muBestTrkType_sublead = 0;
   muTrkNormChi2_sublead = 0;
   muIDPatMVA_sublead = 0;
   muIDPatSoftMVA_sublead = 0;
   muIDSelectorBit_sublead = 0;
   muSvChi2 = 0;
   muSvNDOF = 0;
   muSvX = 0;
   muSvY = 0;
   muSvZ = 0;
   muSvXError = 0;
   muSvYError = 0;
   muSvZError = 0;
   muSvMass = 0;
   muSvCtxy = 0;
   muSvCosAngle = 0;
   kaonMMCharge_lead = 0;
   kaonMMD0_lead = 0;
   kaonMMDz_lead = 0;
   kaonMMD0Error_lead = 0;
   kaonMMDzError_lead = 0;
   kaonMMPt_lead = 0;
   kaonMMEta_lead = 0;
   kaonMMPhi_lead = 0;
   kaonMMVx_lead = 0;
   kaonMMVy_lead = 0;
   kaonMMVz_lead = 0;
   kaonMMEn_lead = 0;
   kaonMMTrkChi2_lead = 0;
   kaonMMTrkNDOF_lead = 0;
   kaonMMTrkNormChi2_lead = 0;
   kaonMMCharge_sublead = 0;
   kaonMMD0_sublead = 0;
   kaonMMDz_sublead = 0;
   kaonMMD0Error_sublead = 0;
   kaonMMDzError_sublead = 0;
   kaonMMPt_sublead = 0;
   kaonMMEta_sublead = 0;
   kaonMMPhi_sublead = 0;
   kaonMMVx_sublead = 0;
   kaonMMVy_sublead = 0;
   kaonMMVz_sublead = 0;
   kaonMMEn_sublead = 0;
   kaonMMTrkChi2_sublead = 0;
   kaonMMTrkNDOF_sublead = 0;
   kaonMMTrkNormChi2_sublead = 0;
   bsMMdRmu = 0;
   bsMMdRkaon = 0;
   bsMMdRJpsiPhi = 0;
   bsMMJpsiMass = 0;
   bsMMPhiMass = 0;
   bsMMBsMass = 0;
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
   //fChain->SetBranchAddress("npfPho", &npfPho, &b_npfPho);
   //fChain->SetBranchAddress("pfphoEt", &pfphoEt, &b_pfphoEt);
   //fChain->SetBranchAddress("pfphoEta", &pfphoEta, &b_pfphoEta);
   //fChain->SetBranchAddress("pfphoPhi", &pfphoPhi, &b_pfphoPhi);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge_lead", &eleCharge_lead, &b_eleCharge_lead);
   fChain->SetBranchAddress("eleChargeConsistent_lead", &eleChargeConsistent_lead, &b_eleChargeConsistent_lead);
   fChain->SetBranchAddress("eleEn_lead", &eleEn_lead, &b_eleEn_lead);
   fChain->SetBranchAddress("eleSCEn_lead", &eleSCEn_lead, &b_eleSCEn_lead);
   fChain->SetBranchAddress("eleEcalEn_lead", &eleEcalEn_lead, &b_eleEcalEn_lead);
   fChain->SetBranchAddress("eleESEnP1_lead", &eleESEnP1_lead, &b_eleESEnP1_lead);
   fChain->SetBranchAddress("eleESEnP2_lead", &eleESEnP2_lead, &b_eleESEnP2_lead);
   fChain->SetBranchAddress("eleD0_lead", &eleD0_lead, &b_eleD0_lead);
   fChain->SetBranchAddress("eleDz_lead", &eleDz_lead, &b_eleDz_lead);
   fChain->SetBranchAddress("eleD0Error_lead", &eleD0Error_lead, &b_eleD0Error_lead);
   fChain->SetBranchAddress("eleDzError_lead", &eleDzError_lead, &b_eleDzError_lead);
   fChain->SetBranchAddress("eleSIP_lead", &eleSIP_lead, &b_eleSIP_lead);
   fChain->SetBranchAddress("elePt_lead", &elePt_lead, &b_elePt_lead);
   fChain->SetBranchAddress("eleEta_lead", &eleEta_lead, &b_eleEta_lead);
   fChain->SetBranchAddress("elePhi_lead", &elePhi_lead, &b_elePhi_lead);
   fChain->SetBranchAddress("eleR9_lead", &eleR9_lead, &b_eleR9_lead);
   fChain->SetBranchAddress("eleCalibPt_lead", &eleCalibPt_lead, &b_eleCalibPt_lead);
   fChain->SetBranchAddress("eleCalibEn_lead", &eleCalibEn_lead, &b_eleCalibEn_lead);
   fChain->SetBranchAddress("eleSCEta_lead", &eleSCEta_lead, &b_eleSCEta_lead);
   fChain->SetBranchAddress("eleSCPhi_lead", &eleSCPhi_lead, &b_eleSCPhi_lead);
   fChain->SetBranchAddress("eleSCRawEn_lead", &eleSCRawEn_lead, &b_eleSCRawEn_lead);
   fChain->SetBranchAddress("eleSCEtaWidth_lead", &eleSCEtaWidth_lead, &b_eleSCEtaWidth_lead);
   fChain->SetBranchAddress("eleSCPhiWidth_lead", &eleSCPhiWidth_lead, &b_eleSCPhiWidth_lead);
   fChain->SetBranchAddress("eleHoverE_lead", &eleHoverE_lead, &b_eleHoverE_lead);
   fChain->SetBranchAddress("eleEoverP_lead", &eleEoverP_lead, &b_eleEoverP_lead);
   fChain->SetBranchAddress("eleEoverPout_lead", &eleEoverPout_lead, &b_eleEoverPout_lead);
   fChain->SetBranchAddress("eleEoverPInv_lead", &eleEoverPInv_lead, &b_eleEoverPInv_lead);
   fChain->SetBranchAddress("eleBrem_lead", &eleBrem_lead, &b_eleBrem_lead);
   fChain->SetBranchAddress("eledEtaAtVtx_lead", &eledEtaAtVtx_lead, &b_eledEtaAtVtx_lead);
   fChain->SetBranchAddress("eledPhiAtVtx_lead", &eledPhiAtVtx_lead, &b_eledPhiAtVtx_lead);
   fChain->SetBranchAddress("eledEtaAtCalo_lead", &eledEtaAtCalo_lead, &b_eledEtaAtCalo_lead);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5_lead", &eleSigmaIEtaIEtaFull5x5_lead, &b_eleSigmaIEtaIEtaFull5x5_lead);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5_lead", &eleSigmaIPhiIPhiFull5x5_lead, &b_eleSigmaIPhiIPhiFull5x5_lead);
   fChain->SetBranchAddress("eleConvVeto_lead", &eleConvVeto_lead, &b_eleConvVeto_lead);
   fChain->SetBranchAddress("eleMissHits_lead", &eleMissHits_lead, &b_eleMissHits_lead);
   fChain->SetBranchAddress("eleESEffSigmaRR_lead", &eleESEffSigmaRR_lead, &b_eleESEffSigmaRR_lead);
   fChain->SetBranchAddress("elePFChIso_lead", &elePFChIso_lead, &b_elePFChIso_lead);
   fChain->SetBranchAddress("elePFPhoIso_lead", &elePFPhoIso_lead, &b_elePFPhoIso_lead);
   fChain->SetBranchAddress("elePFNeuIso_lead", &elePFNeuIso_lead, &b_elePFNeuIso_lead);
   fChain->SetBranchAddress("elePFPUIso_lead", &elePFPUIso_lead, &b_elePFPUIso_lead);
   fChain->SetBranchAddress("elePFClusEcalIso_lead", &elePFClusEcalIso_lead, &b_elePFClusEcalIso_lead);
   fChain->SetBranchAddress("elePFClusHcalIso_lead", &elePFClusHcalIso_lead, &b_elePFClusHcalIso_lead);
   fChain->SetBranchAddress("eleIDMVAIso_lead", &eleIDMVAIso_lead, &b_eleIDMVAIso_lead);
   fChain->SetBranchAddress("eleIDMVANoIso_lead", &eleIDMVANoIso_lead, &b_eleIDMVANoIso_lead);
   fChain->SetBranchAddress("eledEtaseedAtVtx_lead", &eledEtaseedAtVtx_lead, &b_eledEtaseedAtVtx_lead);
   fChain->SetBranchAddress("eleE1x5_lead", &eleE1x5_lead, &b_eleE1x5_lead);
   fChain->SetBranchAddress("eleE2x5_lead", &eleE2x5_lead, &b_eleE2x5_lead);
   fChain->SetBranchAddress("eleE5x5_lead", &eleE5x5_lead, &b_eleE5x5_lead);
   fChain->SetBranchAddress("eleE1x5Full5x5_lead", &eleE1x5Full5x5_lead, &b_eleE1x5Full5x5_lead);
   fChain->SetBranchAddress("eleE2x5Full5x5_lead", &eleE2x5Full5x5_lead, &b_eleE2x5Full5x5_lead);
   fChain->SetBranchAddress("eleE5x5Full5x5_lead", &eleE5x5Full5x5_lead, &b_eleE5x5Full5x5_lead);
   fChain->SetBranchAddress("eleR9Full5x5_lead", &eleR9Full5x5_lead, &b_eleR9Full5x5_lead);
   fChain->SetBranchAddress("eleEcalDrivenSeed_lead", &eleEcalDrivenSeed_lead, &b_eleEcalDrivenSeed_lead);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt_lead", &eleDr03EcalRecHitSumEt_lead, &b_eleDr03EcalRecHitSumEt_lead);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt_lead", &eleDr03HcalDepth1TowerSumEt_lead, &b_eleDr03HcalDepth1TowerSumEt_lead);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt_lead", &eleDr03HcalDepth2TowerSumEt_lead, &b_eleDr03HcalDepth2TowerSumEt_lead);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt_lead", &eleDr03HcalTowerSumEt_lead, &b_eleDr03HcalTowerSumEt_lead);
   fChain->SetBranchAddress("eleDr03TkSumPt_lead", &eleDr03TkSumPt_lead, &b_eleDr03TkSumPt_lead);
   fChain->SetBranchAddress("elecaloEnergy_lead", &elecaloEnergy_lead, &b_elecaloEnergy_lead);
   fChain->SetBranchAddress("eleTrkdxy_lead", &eleTrkdxy_lead, &b_eleTrkdxy_lead);
   fChain->SetBranchAddress("eleKFHits_lead", &eleKFHits_lead, &b_eleKFHits_lead);
   fChain->SetBranchAddress("eleKFChi2_lead", &eleKFChi2_lead, &b_eleKFChi2_lead);
   fChain->SetBranchAddress("eleGSFChi2_lead", &eleGSFChi2_lead, &b_eleGSFChi2_lead);
   fChain->SetBranchAddress("eleIDbit_lead", &eleIDbit_lead, &b_eleIDbit_lead);
   fChain->SetBranchAddress("eleCharge_sublead", &eleCharge_sublead, &b_eleCharge_sublead);
   fChain->SetBranchAddress("eleChargeConsistent_sublead", &eleChargeConsistent_sublead, &b_eleChargeConsistent_sublead);
   fChain->SetBranchAddress("eleEn_sublead", &eleEn_sublead, &b_eleEn_sublead);
   fChain->SetBranchAddress("eleSCEn_sublead", &eleSCEn_sublead, &b_eleSCEn_sublead);
   fChain->SetBranchAddress("eleEcalEn_sublead", &eleEcalEn_sublead, &b_eleEcalEn_sublead);
   fChain->SetBranchAddress("eleESEnP1_sublead", &eleESEnP1_sublead, &b_eleESEnP1_sublead);
   fChain->SetBranchAddress("eleESEnP2_sublead", &eleESEnP2_sublead, &b_eleESEnP2_sublead);
   fChain->SetBranchAddress("eleD0_sublead", &eleD0_sublead, &b_eleD0_sublead);
   fChain->SetBranchAddress("eleDz_sublead", &eleDz_sublead, &b_eleDz_sublead);
   fChain->SetBranchAddress("eleD0Error_sublead", &eleD0Error_sublead, &b_eleD0Error_sublead);
   fChain->SetBranchAddress("eleDzError_sublead", &eleDzError_sublead, &b_eleDzError_sublead);
   fChain->SetBranchAddress("eleSIP_sublead", &eleSIP_sublead, &b_eleSIP_sublead);
   fChain->SetBranchAddress("elePt_sublead", &elePt_sublead, &b_elePt_sublead);
   fChain->SetBranchAddress("eleEta_sublead", &eleEta_sublead, &b_eleEta_sublead);
   fChain->SetBranchAddress("elePhi_sublead", &elePhi_sublead, &b_elePhi_sublead);
   fChain->SetBranchAddress("eleR9_sublead", &eleR9_sublead, &b_eleR9_sublead);
   fChain->SetBranchAddress("eleCalibPt_sublead", &eleCalibPt_sublead, &b_eleCalibPt_sublead);
   fChain->SetBranchAddress("eleCalibEn_sublead", &eleCalibEn_sublead, &b_eleCalibEn_sublead);
   fChain->SetBranchAddress("eleSCEta_sublead", &eleSCEta_sublead, &b_eleSCEta_sublead);
   fChain->SetBranchAddress("eleSCPhi_sublead", &eleSCPhi_sublead, &b_eleSCPhi_sublead);
   fChain->SetBranchAddress("eleSCRawEn_sublead", &eleSCRawEn_sublead, &b_eleSCRawEn_sublead);
   fChain->SetBranchAddress("eleSCEtaWidth_sublead", &eleSCEtaWidth_sublead, &b_eleSCEtaWidth_sublead);
   fChain->SetBranchAddress("eleSCPhiWidth_sublead", &eleSCPhiWidth_sublead, &b_eleSCPhiWidth_sublead);
   fChain->SetBranchAddress("eleHoverE_sublead", &eleHoverE_sublead, &b_eleHoverE_sublead);
   fChain->SetBranchAddress("eleEoverP_sublead", &eleEoverP_sublead, &b_eleEoverP_sublead);
   fChain->SetBranchAddress("eleEoverPout_sublead", &eleEoverPout_sublead, &b_eleEoverPout_sublead);
   fChain->SetBranchAddress("eleEoverPInv_sublead", &eleEoverPInv_sublead, &b_eleEoverPInv_sublead);
   fChain->SetBranchAddress("eleBrem_sublead", &eleBrem_sublead, &b_eleBrem_sublead);
   fChain->SetBranchAddress("eledEtaAtVtx_sublead", &eledEtaAtVtx_sublead, &b_eledEtaAtVtx_sublead);
   fChain->SetBranchAddress("eledPhiAtVtx_sublead", &eledPhiAtVtx_sublead, &b_eledPhiAtVtx_sublead);
   fChain->SetBranchAddress("eledEtaAtCalo_sublead", &eledEtaAtCalo_sublead, &b_eledEtaAtCalo_sublead);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5_sublead", &eleSigmaIEtaIEtaFull5x5_sublead, &b_eleSigmaIEtaIEtaFull5x5_sublead);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5_sublead", &eleSigmaIPhiIPhiFull5x5_sublead, &b_eleSigmaIPhiIPhiFull5x5_sublead);
   fChain->SetBranchAddress("eleConvVeto_sublead", &eleConvVeto_sublead, &b_eleConvVeto_sublead);
   fChain->SetBranchAddress("eleMissHits_sublead", &eleMissHits_sublead, &b_eleMissHits_sublead);
   fChain->SetBranchAddress("eleESEffSigmaRR_sublead", &eleESEffSigmaRR_sublead, &b_eleESEffSigmaRR_sublead);
   fChain->SetBranchAddress("elePFChIso_sublead", &elePFChIso_sublead, &b_elePFChIso_sublead);
   fChain->SetBranchAddress("elePFPhoIso_sublead", &elePFPhoIso_sublead, &b_elePFPhoIso_sublead);
   fChain->SetBranchAddress("elePFNeuIso_sublead", &elePFNeuIso_sublead, &b_elePFNeuIso_sublead);
   fChain->SetBranchAddress("elePFPUIso_sublead", &elePFPUIso_sublead, &b_elePFPUIso_sublead);
   fChain->SetBranchAddress("elePFClusEcalIso_sublead", &elePFClusEcalIso_sublead, &b_elePFClusEcalIso_sublead);
   fChain->SetBranchAddress("elePFClusHcalIso_sublead", &elePFClusHcalIso_sublead, &b_elePFClusHcalIso_sublead);
   fChain->SetBranchAddress("eleIDMVAIso_sublead", &eleIDMVAIso_sublead, &b_eleIDMVAIso_sublead);
   fChain->SetBranchAddress("eleIDMVANoIso_sublead", &eleIDMVANoIso_sublead, &b_eleIDMVANoIso_sublead);
   fChain->SetBranchAddress("eledEtaseedAtVtx_sublead", &eledEtaseedAtVtx_sublead, &b_eledEtaseedAtVtx_sublead);
   fChain->SetBranchAddress("eleE1x5_sublead", &eleE1x5_sublead, &b_eleE1x5_sublead);
   fChain->SetBranchAddress("eleE2x5_sublead", &eleE2x5_sublead, &b_eleE2x5_sublead);
   fChain->SetBranchAddress("eleE5x5_sublead", &eleE5x5_sublead, &b_eleE5x5_sublead);
   fChain->SetBranchAddress("eleE1x5Full5x5_sublead", &eleE1x5Full5x5_sublead, &b_eleE1x5Full5x5_sublead);
   fChain->SetBranchAddress("eleE2x5Full5x5_sublead", &eleE2x5Full5x5_sublead, &b_eleE2x5Full5x5_sublead);
   fChain->SetBranchAddress("eleE5x5Full5x5_sublead", &eleE5x5Full5x5_sublead, &b_eleE5x5Full5x5_sublead);
   fChain->SetBranchAddress("eleR9Full5x5_sublead", &eleR9Full5x5_sublead, &b_eleR9Full5x5_sublead);
   fChain->SetBranchAddress("eleEcalDrivenSeed_sublead", &eleEcalDrivenSeed_sublead, &b_eleEcalDrivenSeed_sublead);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt_sublead", &eleDr03EcalRecHitSumEt_sublead, &b_eleDr03EcalRecHitSumEt_sublead);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt_sublead", &eleDr03HcalDepth1TowerSumEt_sublead, &b_eleDr03HcalDepth1TowerSumEt_sublead);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt_sublead", &eleDr03HcalDepth2TowerSumEt_sublead, &b_eleDr03HcalDepth2TowerSumEt_sublead);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt_sublead", &eleDr03HcalTowerSumEt_sublead, &b_eleDr03HcalTowerSumEt_sublead);
   fChain->SetBranchAddress("eleDr03TkSumPt_sublead", &eleDr03TkSumPt_sublead, &b_eleDr03TkSumPt_sublead);
   fChain->SetBranchAddress("elecaloEnergy_sublead", &elecaloEnergy_sublead, &b_elecaloEnergy_sublead);
   fChain->SetBranchAddress("eleTrkdxy_sublead", &eleTrkdxy_sublead, &b_eleTrkdxy_sublead);
   fChain->SetBranchAddress("eleKFHits_sublead", &eleKFHits_sublead, &b_eleKFHits_sublead);
   fChain->SetBranchAddress("eleKFChi2_sublead", &eleKFChi2_sublead, &b_eleKFChi2_sublead);
   fChain->SetBranchAddress("eleGSFChi2_sublead", &eleGSFChi2_sublead, &b_eleGSFChi2_sublead);
   fChain->SetBranchAddress("eleIDbit_sublead", &eleIDbit_sublead, &b_eleIDbit_sublead);
   fChain->SetBranchAddress("eleSvChi2", &eleSvChi2, &b_eleSvChi2);
   fChain->SetBranchAddress("eleSvNDOF", &eleSvNDOF, &b_eleSvNDOF);
   fChain->SetBranchAddress("eleSvProb", &eleSvProb, &b_eleSvProb);
   fChain->SetBranchAddress("eleSvX", &eleSvX, &b_eleSvX);
   fChain->SetBranchAddress("eleSvY", &eleSvY, &b_eleSvY);
   fChain->SetBranchAddress("eleSvZ", &eleSvZ, &b_eleSvZ);
   fChain->SetBranchAddress("eleSvXError", &eleSvXError, &b_eleSvXError);
   fChain->SetBranchAddress("eleSvYError", &eleSvYError, &b_eleSvYError);
   fChain->SetBranchAddress("eleSvZError", &eleSvZError, &b_eleSvZError);
   fChain->SetBranchAddress("eleSvMass", &eleSvMass, &b_eleSvMass);
   fChain->SetBranchAddress("eleSvCtxy", &eleSvCtxy, &b_eleSvCtxy);
   fChain->SetBranchAddress("eleSvCosAngle", &eleSvCosAngle, &b_eleSvCosAngle);
   fChain->SetBranchAddress("eleSvLxy", &eleSvLxy, &b_eleSvLxy);
   fChain->SetBranchAddress("eleSvLxyError", &eleSvLxyError, &b_eleSvLxyError);
   fChain->SetBranchAddress("kaonEECharge_lead", &kaonEECharge_lead, &b_kaonEECharge_lead);
   fChain->SetBranchAddress("kaonEED0_lead", &kaonEED0_lead, &b_kaonEED0_lead);
   fChain->SetBranchAddress("kaonEEDz_lead", &kaonEEDz_lead, &b_kaonEEDz_lead);
   fChain->SetBranchAddress("kaonEED0Error_lead", &kaonEED0Error_lead, &b_kaonEED0Error_lead);
   fChain->SetBranchAddress("kaonEEDzError_lead", &kaonEEDzError_lead, &b_kaonEEDzError_lead);
   fChain->SetBranchAddress("kaonEEPt_lead", &kaonEEPt_lead, &b_kaonEEPt_lead);
   fChain->SetBranchAddress("kaonEEEta_lead", &kaonEEEta_lead, &b_kaonEEEta_lead);
   fChain->SetBranchAddress("kaonEEPhi_lead", &kaonEEPhi_lead, &b_kaonEEPhi_lead);
   fChain->SetBranchAddress("kaonEEVx_lead", &kaonEEVx_lead, &b_kaonEEVx_lead);
   fChain->SetBranchAddress("kaonEEVy_lead", &kaonEEVy_lead, &b_kaonEEVy_lead);
   fChain->SetBranchAddress("kaonEEVz_lead", &kaonEEVz_lead, &b_kaonEEVz_lead);
   fChain->SetBranchAddress("kaonEETrkChi2_lead", &kaonEETrkChi2_lead, &b_kaonEETrkChi2_lead);
   fChain->SetBranchAddress("kaonEETrkNDOF_lead", &kaonEETrkNDOF_lead, &b_kaonEETrkNDOF_lead);
   fChain->SetBranchAddress("kaonEETrkNormChi2_lead", &kaonEETrkNormChi2_lead, &b_kaonEETrkNormChi2_lead);
   fChain->SetBranchAddress("kaonEECharge_sublead", &kaonEECharge_sublead, &b_kaonEECharge_sublead);
   fChain->SetBranchAddress("kaonEED0_sublead", &kaonEED0_sublead, &b_kaonEED0_sublead);
   fChain->SetBranchAddress("kaonEEDz_sublead", &kaonEEDz_sublead, &b_kaonEEDz_sublead);
   fChain->SetBranchAddress("kaonEED0Error_sublead", &kaonEED0Error_sublead, &b_kaonEED0Error_sublead);
   fChain->SetBranchAddress("kaonEEDzError_sublead", &kaonEEDzError_sublead, &b_kaonEEDzError_sublead);
   fChain->SetBranchAddress("kaonEEPt_sublead", &kaonEEPt_sublead, &b_kaonEEPt_sublead);
   fChain->SetBranchAddress("kaonEEEta_sublead", &kaonEEEta_sublead, &b_kaonEEEta_sublead);
   fChain->SetBranchAddress("kaonEEPhi_sublead", &kaonEEPhi_sublead, &b_kaonEEPhi_sublead);
   fChain->SetBranchAddress("kaonEEVx_sublead", &kaonEEVx_sublead, &b_kaonEEVx_sublead);
   fChain->SetBranchAddress("kaonEEVy_sublead", &kaonEEVy_sublead, &b_kaonEEVy_sublead);
   fChain->SetBranchAddress("kaonEEVz_sublead", &kaonEEVz_sublead, &b_kaonEEVz_sublead);
   fChain->SetBranchAddress("kaonEETrkChi2_sublead", &kaonEETrkChi2_sublead, &b_kaonEETrkChi2_sublead);
   fChain->SetBranchAddress("kaonEETrkNDOF_sublead", &kaonEETrkNDOF_sublead, &b_kaonEETrkNDOF_sublead);
   fChain->SetBranchAddress("kaonEETrkNormChi2_sublead", &kaonEETrkNormChi2_sublead, &b_kaonEETrkNormChi2_sublead);
   fChain->SetBranchAddress("bsEEdRele", &bsEEdRele, &b_bsEEdRele);
   fChain->SetBranchAddress("bsEEdRkaon", &bsEEdRkaon, &b_bsEEdRkaon);
   fChain->SetBranchAddress("bsEEdRJpsiPhi", &bsEEdRJpsiPhi, &b_bsEEdRJpsiPhi);
   fChain->SetBranchAddress("bsEEJpsiMass", &bsEEJpsiMass, &b_bsEEJpsiMass);
   fChain->SetBranchAddress("bsEEPhiMass", &bsEEPhiMass, &b_bsEEPhiMass);
   fChain->SetBranchAddress("bsEEBsMass", &bsEEBsMass, &b_bsEEBsMass);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("matchedTrg", &matchedTrg, &b_matchedTrg);
   fChain->SetBranchAddress("muPt_lead", &muPt_lead, &b_muPt_lead);
   fChain->SetBranchAddress("muEn_lead", &muEn_lead, &b_muEn_lead);
   fChain->SetBranchAddress("muEta_lead", &muEta_lead, &b_muEta_lead);
   fChain->SetBranchAddress("muPhi_lead", &muPhi_lead, &b_muPhi_lead);
   fChain->SetBranchAddress("muCharge_lead", &muCharge_lead, &b_muCharge_lead);
   fChain->SetBranchAddress("muType_lead", &muType_lead, &b_muType_lead);
   fChain->SetBranchAddress("muIDbit_lead", &muIDbit_lead, &b_muIDbit_lead);
   fChain->SetBranchAddress("muD0_lead", &muD0_lead, &b_muD0_lead);
   fChain->SetBranchAddress("muDz_lead", &muDz_lead, &b_muDz_lead);
   fChain->SetBranchAddress("muSIP_lead", &muSIP_lead, &b_muSIP_lead);
   fChain->SetBranchAddress("muD0Error_lead", &muD0Error_lead, &b_muD0Error_lead);
   fChain->SetBranchAddress("muDzError_lead", &muDzError_lead, &b_muDzError_lead);
   fChain->SetBranchAddress("muChi2NDF_lead", &muChi2NDF_lead, &b_muChi2NDF_lead);
   fChain->SetBranchAddress("muInnerD0_lead", &muInnerD0_lead, &b_muInnerD0_lead);
   fChain->SetBranchAddress("muInnerDz_lead", &muInnerDz_lead, &b_muInnerDz_lead);
   fChain->SetBranchAddress("muTrkLayers_lead", &muTrkLayers_lead, &b_muTrkLayers_lead);
   fChain->SetBranchAddress("muPixelLayers_lead", &muPixelLayers_lead, &b_muPixelLayers_lead);
   fChain->SetBranchAddress("muPixelHits_lead", &muPixelHits_lead, &b_muPixelHits_lead);
   fChain->SetBranchAddress("muMuonHits_lead", &muMuonHits_lead, &b_muMuonHits_lead);
   fChain->SetBranchAddress("muStations_lead", &muStations_lead, &b_muStations_lead);
   fChain->SetBranchAddress("muMatches_lead", &muMatches_lead, &b_muMatches_lead);
   fChain->SetBranchAddress("muTrkQuality_lead", &muTrkQuality_lead, &b_muTrkQuality_lead);
   fChain->SetBranchAddress("muIsoTrk_lead", &muIsoTrk_lead, &b_muIsoTrk_lead);
   fChain->SetBranchAddress("muPFChIso_lead", &muPFChIso_lead, &b_muPFChIso_lead);
   fChain->SetBranchAddress("muPFPhoIso_lead", &muPFPhoIso_lead, &b_muPFPhoIso_lead);
   fChain->SetBranchAddress("muPFNeuIso_lead", &muPFNeuIso_lead, &b_muPFNeuIso_lead);
   fChain->SetBranchAddress("muPFPUIso_lead", &muPFPUIso_lead, &b_muPFPUIso_lead);
   fChain->SetBranchAddress("muPFChIso03_lead", &muPFChIso03_lead, &b_muPFChIso03_lead);
   fChain->SetBranchAddress("muPFPhoIso03_lead", &muPFPhoIso03_lead, &b_muPFPhoIso03_lead);
   fChain->SetBranchAddress("muPFNeuIso03_lead", &muPFNeuIso03_lead, &b_muPFNeuIso03_lead);
   fChain->SetBranchAddress("muPFPUIso03_lead", &muPFPUIso03_lead, &b_muPFPUIso03_lead);
   fChain->SetBranchAddress("muFiredTrgs_lead", &muFiredTrgs_lead, &b_muFiredTrgs_lead);
   fChain->SetBranchAddress("muFiredL1Trgs_lead", &muFiredL1Trgs_lead, &b_muFiredL1Trgs_lead);
   fChain->SetBranchAddress("muInnervalidFraction_lead", &muInnervalidFraction_lead, &b_muInnervalidFraction_lead);
   fChain->SetBranchAddress("musegmentCompatibility_lead", &musegmentCompatibility_lead, &b_musegmentCompatibility_lead);
   fChain->SetBranchAddress("muchi2LocalPosition_lead", &muchi2LocalPosition_lead, &b_muchi2LocalPosition_lead);
   fChain->SetBranchAddress("mutrkKink_lead", &mutrkKink_lead, &b_mutrkKink_lead);
   fChain->SetBranchAddress("muBestTrkPtError_lead", &muBestTrkPtError_lead, &b_muBestTrkPtError_lead);
   fChain->SetBranchAddress("muBestTrkPt_lead", &muBestTrkPt_lead, &b_muBestTrkPt_lead);
   fChain->SetBranchAddress("muBestTrkType_lead", &muBestTrkType_lead, &b_muBestTrkType_lead);
   fChain->SetBranchAddress("muTrkNormChi2_lead", &muTrkNormChi2_lead, &b_muTrkNormChi2_lead);
   fChain->SetBranchAddress("muIDPatMVA_lead", &muIDPatMVA_lead, &b_muIDPatMVA_lead);
   fChain->SetBranchAddress("muIDPatSoftMVA_lead", &muIDPatSoftMVA_lead, &b_muIDPatSoftMVA_lead);
   fChain->SetBranchAddress("muIDSelectorBit_lead", &muIDSelectorBit_lead, &b_muIDSelectorBit_lead);
   fChain->SetBranchAddress("muPt_sublead", &muPt_sublead, &b_muPt_sublead);
   fChain->SetBranchAddress("muEn_sublead", &muEn_sublead, &b_muEn_sublead);
   fChain->SetBranchAddress("muEta_sublead", &muEta_sublead, &b_muEta_sublead);
   fChain->SetBranchAddress("muPhi_sublead", &muPhi_sublead, &b_muPhi_sublead);
   fChain->SetBranchAddress("muCharge_sublead", &muCharge_sublead, &b_muCharge_sublead);
   fChain->SetBranchAddress("muType_sublead", &muType_sublead, &b_muType_sublead);
   fChain->SetBranchAddress("muIDbit_sublead", &muIDbit_sublead, &b_muIDbit_sublead);
   fChain->SetBranchAddress("muD0_sublead", &muD0_sublead, &b_muD0_sublead);
   fChain->SetBranchAddress("muDz_sublead", &muDz_sublead, &b_muDz_sublead);
   fChain->SetBranchAddress("muSIP_sublead", &muSIP_sublead, &b_muSIP_sublead);
   fChain->SetBranchAddress("muD0Error_sublead", &muD0Error_sublead, &b_muD0Error_sublead);
   fChain->SetBranchAddress("muDzError_sublead", &muDzError_sublead, &b_muDzError_sublead);
   fChain->SetBranchAddress("muChi2NDF_sublead", &muChi2NDF_sublead, &b_muChi2NDF_sublead);
   fChain->SetBranchAddress("muInnerD0_sublead", &muInnerD0_sublead, &b_muInnerD0_sublead);
   fChain->SetBranchAddress("muInnerDz_sublead", &muInnerDz_sublead, &b_muInnerDz_sublead);
   fChain->SetBranchAddress("muTrkLayers_sublead", &muTrkLayers_sublead, &b_muTrkLayers_sublead);
   fChain->SetBranchAddress("muPixelLayers_sublead", &muPixelLayers_sublead, &b_muPixelLayers_sublead);
   fChain->SetBranchAddress("muPixelHits_sublead", &muPixelHits_sublead, &b_muPixelHits_sublead);
   fChain->SetBranchAddress("muMuonHits_sublead", &muMuonHits_sublead, &b_muMuonHits_sublead);
   fChain->SetBranchAddress("muStations_sublead", &muStations_sublead, &b_muStations_sublead);
   fChain->SetBranchAddress("muMatches_sublead", &muMatches_sublead, &b_muMatches_sublead);
   fChain->SetBranchAddress("muTrkQuality_sublead", &muTrkQuality_sublead, &b_muTrkQuality_sublead);
   fChain->SetBranchAddress("muIsoTrk_sublead", &muIsoTrk_sublead, &b_muIsoTrk_sublead);
   fChain->SetBranchAddress("muPFChIso_sublead", &muPFChIso_sublead, &b_muPFChIso_sublead);
   fChain->SetBranchAddress("muPFPhoIso_sublead", &muPFPhoIso_sublead, &b_muPFPhoIso_sublead);
   fChain->SetBranchAddress("muPFNeuIso_sublead", &muPFNeuIso_sublead, &b_muPFNeuIso_sublead);
   fChain->SetBranchAddress("muPFPUIso_sublead", &muPFPUIso_sublead, &b_muPFPUIso_sublead);
   fChain->SetBranchAddress("muPFChIso03_sublead", &muPFChIso03_sublead, &b_muPFChIso03_sublead);
   fChain->SetBranchAddress("muPFPhoIso03_sublead", &muPFPhoIso03_sublead, &b_muPFPhoIso03_sublead);
   fChain->SetBranchAddress("muPFNeuIso03_sublead", &muPFNeuIso03_sublead, &b_muPFNeuIso03_sublead);
   fChain->SetBranchAddress("muPFPUIso03_sublead", &muPFPUIso03_sublead, &b_muPFPUIso03_sublead);
   fChain->SetBranchAddress("muFiredTrgs_sublead", &muFiredTrgs_sublead, &b_muFiredTrgs_sublead);
   fChain->SetBranchAddress("muFiredL1Trgs_sublead", &muFiredL1Trgs_sublead, &b_muFiredL1Trgs_sublead);
   fChain->SetBranchAddress("muInnervalidFraction_sublead", &muInnervalidFraction_sublead, &b_muInnervalidFraction_sublead);
   fChain->SetBranchAddress("musegmentCompatibility_sublead", &musegmentCompatibility_sublead, &b_musegmentCompatibility_sublead);
   fChain->SetBranchAddress("muchi2LocalPosition_sublead", &muchi2LocalPosition_sublead, &b_muchi2LocalPosition_sublead);
   fChain->SetBranchAddress("mutrkKink_sublead", &mutrkKink_sublead, &b_mutrkKink_sublead);
   fChain->SetBranchAddress("muBestTrkPtError_sublead", &muBestTrkPtError_sublead, &b_muBestTrkPtError_sublead);
   fChain->SetBranchAddress("muBestTrkPt_sublead", &muBestTrkPt_sublead, &b_muBestTrkPt_sublead);
   fChain->SetBranchAddress("muBestTrkType_sublead", &muBestTrkType_sublead, &b_muBestTrkType_sublead);
   fChain->SetBranchAddress("muTrkNormChi2_sublead", &muTrkNormChi2_sublead, &b_muTrkNormChi2_sublead);
   fChain->SetBranchAddress("muIDPatMVA_sublead", &muIDPatMVA_sublead, &b_muIDPatMVA_sublead);
   fChain->SetBranchAddress("muIDPatSoftMVA_sublead", &muIDPatSoftMVA_sublead, &b_muIDPatSoftMVA_sublead);
   fChain->SetBranchAddress("muIDSelectorBit_sublead", &muIDSelectorBit_sublead, &b_muIDSelectorBit_sublead);
   fChain->SetBranchAddress("muSvChi2", &muSvChi2, &b_muSvChi2);
   fChain->SetBranchAddress("muSvNDOF", &muSvNDOF, &b_muSvNDOF);
   fChain->SetBranchAddress("muSvX", &muSvX, &b_muSvX);
   fChain->SetBranchAddress("muSvY", &muSvY, &b_muSvY);
   fChain->SetBranchAddress("muSvZ", &muSvZ, &b_muSvZ);
   fChain->SetBranchAddress("muSvXError", &muSvXError, &b_muSvXError);
   fChain->SetBranchAddress("muSvYError", &muSvYError, &b_muSvYError);
   fChain->SetBranchAddress("muSvZError", &muSvZError, &b_muSvZError);
   fChain->SetBranchAddress("muSvMass", &muSvMass, &b_muSvMass);
   fChain->SetBranchAddress("muSvCtxy", &muSvCtxy, &b_muSvCtxy);
   fChain->SetBranchAddress("muSvCosAngle", &muSvCosAngle, &b_muSvCosAngle);
   fChain->SetBranchAddress("kaonMMCharge_lead", &kaonMMCharge_lead, &b_kaonMMCharge_lead);
   fChain->SetBranchAddress("kaonMMD0_lead", &kaonMMD0_lead, &b_kaonMMD0_lead);
   fChain->SetBranchAddress("kaonMMDz_lead", &kaonMMDz_lead, &b_kaonMMDz_lead);
   fChain->SetBranchAddress("kaonMMD0Error_lead", &kaonMMD0Error_lead, &b_kaonMMD0Error_lead);
   fChain->SetBranchAddress("kaonMMDzError_lead", &kaonMMDzError_lead, &b_kaonMMDzError_lead);
   fChain->SetBranchAddress("kaonMMPt_lead", &kaonMMPt_lead, &b_kaonMMPt_lead);
   fChain->SetBranchAddress("kaonMMEta_lead", &kaonMMEta_lead, &b_kaonMMEta_lead);
   fChain->SetBranchAddress("kaonMMPhi_lead", &kaonMMPhi_lead, &b_kaonMMPhi_lead);
   fChain->SetBranchAddress("kaonMMVx_lead", &kaonMMVx_lead, &b_kaonMMVx_lead);
   fChain->SetBranchAddress("kaonMMVy_lead", &kaonMMVy_lead, &b_kaonMMVy_lead);
   fChain->SetBranchAddress("kaonMMVz_lead", &kaonMMVz_lead, &b_kaonMMVz_lead);
   fChain->SetBranchAddress("kaonMMEn_lead", &kaonMMEn_lead, &b_kaonMMEn_lead);
   fChain->SetBranchAddress("kaonMMTrkChi2_lead", &kaonMMTrkChi2_lead, &b_kaonMMTrkChi2_lead);
   fChain->SetBranchAddress("kaonMMTrkNDOF_lead", &kaonMMTrkNDOF_lead, &b_kaonMMTrkNDOF_lead);
   fChain->SetBranchAddress("kaonMMTrkNormChi2_lead", &kaonMMTrkNormChi2_lead, &b_kaonMMTrkNormChi2_lead);
   fChain->SetBranchAddress("kaonMMCharge_sublead", &kaonMMCharge_sublead, &b_kaonMMCharge_sublead);
   fChain->SetBranchAddress("kaonMMD0_sublead", &kaonMMD0_sublead, &b_kaonMMD0_sublead);
   fChain->SetBranchAddress("kaonMMDz_sublead", &kaonMMDz_sublead, &b_kaonMMDz_sublead);
   fChain->SetBranchAddress("kaonMMD0Error_sublead", &kaonMMD0Error_sublead, &b_kaonMMD0Error_sublead);
   fChain->SetBranchAddress("kaonMMDzError_sublead", &kaonMMDzError_sublead, &b_kaonMMDzError_sublead);
   fChain->SetBranchAddress("kaonMMPt_sublead", &kaonMMPt_sublead, &b_kaonMMPt_sublead);
   fChain->SetBranchAddress("kaonMMEta_sublead", &kaonMMEta_sublead, &b_kaonMMEta_sublead);
   fChain->SetBranchAddress("kaonMMPhi_sublead", &kaonMMPhi_sublead, &b_kaonMMPhi_sublead);
   fChain->SetBranchAddress("kaonMMVx_sublead", &kaonMMVx_sublead, &b_kaonMMVx_sublead);
   fChain->SetBranchAddress("kaonMMVy_sublead", &kaonMMVy_sublead, &b_kaonMMVy_sublead);
   fChain->SetBranchAddress("kaonMMVz_sublead", &kaonMMVz_sublead, &b_kaonMMVz_sublead);
   fChain->SetBranchAddress("kaonMMEn_sublead", &kaonMMEn_sublead, &b_kaonMMEn_sublead);
   fChain->SetBranchAddress("kaonMMTrkChi2_sublead", &kaonMMTrkChi2_sublead, &b_kaonMMTrkChi2_sublead);
   fChain->SetBranchAddress("kaonMMTrkNDOF_sublead", &kaonMMTrkNDOF_sublead, &b_kaonMMTrkNDOF_sublead);
   fChain->SetBranchAddress("kaonMMTrkNormChi2_sublead", &kaonMMTrkNormChi2_sublead, &b_kaonMMTrkNormChi2_sublead);
   fChain->SetBranchAddress("bsMMdRmu", &bsMMdRmu, &b_bsMMdRmu);
   fChain->SetBranchAddress("bsMMdRkaon", &bsMMdRkaon, &b_bsMMdRkaon);
   fChain->SetBranchAddress("bsMMdRJpsiPhi", &bsMMdRJpsiPhi, &b_bsMMdRJpsiPhi);
   fChain->SetBranchAddress("bsMMJpsiMass", &bsMMJpsiMass, &b_bsMMJpsiMass);
   fChain->SetBranchAddress("bsMMPhiMass", &bsMMPhiMass, &b_bsMMPhiMass);
   fChain->SetBranchAddress("bsMMBsMass", &bsMMBsMass, &b_bsMMBsMass);
   Notify();
}

Bool_t BsPhiLLTupleTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BsPhiLLTupleTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BsPhiLLTupleTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Float_t BsPhiLLTupleTree::deltaPhi(Float_t phi1, Float_t phi2)
{
   Float_t dphi = phi1 - phi2;
   while (dphi >= TMath::Pi()) dphi -= 2*TMath::Pi();
   while (dphi < -1*TMath::Pi()) dphi += 2*TMath::Pi();
   return dphi;
}
#endif // #ifdef BsPhiLLTupleTree_cxx
