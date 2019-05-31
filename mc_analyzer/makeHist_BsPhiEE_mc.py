#! /usr/bin/env python

import argparse

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="DoubleMuonNtu_Run2016B.list", help="List of input ggNtuplizer files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-m", "--maxevents", dest="maxevents", type=int, default=-1, help="Maximum number events to loop over")
parser.add_argument("-t", "--ttree", dest="ttree", default="tree", help="TTree Name")
args = parser.parse_args()

import numpy as np
import itertools
import ROOT
import os
from collections import OrderedDict
from array import array


if os.path.isfile('~/.rootlogon.C'): ROOT.gROOT.Macro(os.path.expanduser('~/.rootlogon.C'))
ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(000000)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.UseCurrentStyle()

sw = ROOT.TStopwatch()
sw.Start()

# Input ggNtuple
tchain = ROOT.TChain(args.ttree)
"""
with open(args.inputfiles) as filenames:
    for filename in filenames: 
        tchain.Add(filename.rstrip('\n'))
"""
tchain.Add(args.inputfiles)
print('Total number of events: ' + str(tchain.GetEntries()))

# Output file and any histograms we want
file_out = ROOT.TFile(args.outputfile, 'recreate')

h_q2_InvM_mc = ROOT.TH1D("h_q2_InvM_mc", "", 200, 0.0, 25.0);
h_jpsi_InvM_mc = ROOT.TH1D("h_jpsi_InvM_mc", "", 100, 2.4, 3.8);
h_phiee_InvM_mc = ROOT.TH1D("h_phiee_InvM_mc", "", 200, 0.98, 1.1);
h_bsee_InvM_mc = ROOT.TH1D("h_bsee_InvM_mc", "", 200, 4.5, 6.0);

h_electron_pt_lead_mc = ROOT.TH1D("h_electron_pt_lead_mc", "", 80, 0.0, 20.0);
h_electron_eta_lead_mc = ROOT.TH1D("h_electron_eta_lead_mc", "", 50, -3, 3);
h_electron_phi_lead_mc = ROOT.TH1D("h_electron_phi_lead_mc", "", 50, -1.0*ROOT.TMath.Pi(), ROOT.TMath.Pi());
h_electron_D0_lead_mc = ROOT.TH1D("h_electron_D0_lead_mc", "", 100, -0.5, 0.5);
h_electron_Dz_lead_mc = ROOT.TH1D("h_electron_Dz_lead_mc", "", 100, -1.0, 1.0);

h_electron_pt_sublead_mc = ROOT.TH1D("h_electron_pt_sublead_mc", "", 80, 0.0, 20.0);
h_electron_eta_sublead_mc = ROOT.TH1D("h_electron_eta_sublead_mc", "", 50, -3, 3);
h_electron_phi_sublead_mc = ROOT.TH1D("h_electron_phi_sublead_mc", "", 50, -1.0*ROOT.TMath.Pi(), ROOT.TMath.Pi());
h_electron_D0_sublead_mc = ROOT.TH1D("h_electron_D0_sublead_mc", "", 100, -0.5, 0.5);
h_electron_Dz_sublead_mc = ROOT.TH1D("h_electron_Dz_sublead_mc", "", 100, -1.0, 1.0);

h_kaon_ee_pt_lead_mc = ROOT.TH1D("h_kaon_ee_pt_lead_mc", "", 80, 0.0, 5.0);
h_kaon_ee_eta_lead_mc = ROOT.TH1D("h_kaon_ee_eta_lead_mc", "", 50, -3, 3);
h_kaon_ee_phi_lead_mc = ROOT.TH1D("h_kaon_ee_phi_lead_mc", "", 50, -1.0*ROOT.TMath.Pi(), ROOT.TMath.Pi());
h_kaon_ee_Dz_lead_mc = ROOT.TH1D("h_kaon_ee_Dz_lead_mc", "", 100, -1.0, 1.0);
h_kaon_ee_D0_lead_mc = ROOT.TH1D("h_kaon_ee_D0_lead_mc", "", 100, -0.5, 0.5);
h_kaon_ee_D0Sig_lead_mc = ROOT.TH1D("h_kaon_ee_D0Sig_lead_mc", "", 100, -20.0, 20.0);
h_kaon_ee_normChi2_lead_mc = ROOT.TH1D("h_kaon_ee_normChi2_lead_mc", "", 100, 0.0, 5.0);

h_kaon_ee_pt_sublead_mc = ROOT.TH1D("h_kaon_ee_pt_sublead_mc", "", 80, 0.0, 5.0);
h_kaon_ee_eta_sublead_mc = ROOT.TH1D("h_kaon_ee_eta_sublead_mc", "", 50, -3, 3);
h_kaon_ee_phi_sublead_mc = ROOT.TH1D("h_kaon_ee_phi_sublead_mc", "", 50, -1.0*ROOT.TMath.Pi(), ROOT.TMath.Pi());
h_kaon_ee_Dz_sublead_mc = ROOT.TH1D("h_kaon_ee_Dz_sublead_mc", "", 100, -1.0, 1.0);
h_kaon_ee_D0_sublead_mc = ROOT.TH1D("h_kaon_ee_D0_sublead_mc", "", 100, -0.5, 0.5);
h_kaon_ee_D0Sig_sublead_mc = ROOT.TH1D("h_kaon_ee_D0Sig_sublead_mc", "", 100, -20.0, 20.0);
h_kaon_ee_normChi2_sublead_mc = ROOT.TH1D("h_kaon_ee_normChi2_sublead_mc", "", 100, 0.0, 5.0);

h_bs_pt_ee_mc = ROOT.TH1D("h_bs_pt_ee_mc", "", 100, 0.0, 50.0);

h_svProb_ee_mc = ROOT.TH1D("h_svProb_ee_mc", "", 100, 0.0, 1.0);
h_svCosAngle_ee_mc = ROOT.TH1D("h_svCosAngle_ee_mc", "", 200, -1.0, 1.0);
h_svCtxy_ee_mc = ROOT.TH1D("h_svCtxy_ee_mc", "", 100, -1.0, 1.0);
h_svLxy_ee_mc = ROOT.TH1D("h_svLxy_ee_mc", "", 100, 0.0, 1);
h_svLxySig_ee_mc = ROOT.TH1D("h_svLxySig_ee_mc", "", 200, 0.0, 60);

h_electron_dR_mc = ROOT.TH1D("h_electron_dR_mc", "", 120, 0.0, 4.0);
h_kaon_ee_dR_mc = ROOT.TH1D("h_kaon_ee_dR_mc", "", 120, 0.0, 4.0);
h_jpsiPhiOpen_ee_mc = ROOT.TH1D("h_jpsiPhiOpen_ee_mc", "", 120, 0.0, 4.0);



eleM = 0.0005109989461
kaonM = 0.493677
pionM = 0.13957018
bsM = 5.3663
lv_ele_lead_gen = ROOT.TLorentzVector()
lv_ele_sublead_gen = ROOT.TLorentzVector()
lv_kaon_lead_gen = ROOT.TLorentzVector()
lv_kaon_sublead_gen = ROOT.TLorentzVector()
lv_ele_lead_reco = ROOT.TLorentzVector()
lv_ele_sublead_reco = ROOT.TLorentzVector()
lv_kaon_lead_reco = ROOT.TLorentzVector()
lv_kaon_sublead_reco = ROOT.TLorentzVector()
lv_bs_reco = ROOT.TLorentzVector()

nTot = 0

#Loop over all the events in the input ntuple
for ievent,event in enumerate(tchain):
    if ievent > args.maxevents and args.maxevents != -1: break
    if ievent % 10000 == 0: print('Processing entry ' + str(ievent))

    if event.tag_mu_reco_dxy == -99: continue
    if abs(event.b0_lp_eta) > 2.5 or abs(event.b0_lm_eta) > 2.5 or abs(event.b0_k_eta) > 2.5 or abs(event.b0_pi_eta) > 2.5: continue
    if abs(event.b0_lp_reco_pt) < 0.4 or abs(event.b0_lm_reco_pt) < 0.4 or abs(event.b0_k_reco_pt) < 0.4 or abs(event.b0_pi_reco_pt) < 0.4: continue
    if event.foundB0 < 0.7: continue
    #if event.tag_mu_reco_pt < 9: continue
    #if abs(event.tag_mu_reco_dxy/event.tag_mu_reco_dxyError) < 6: continue

    if event.b0_lp_pt > event.b0_lm_pt:
        lv_ele_lead_reco.SetPtEtaPhiM(event.b0_lp_reco_pt, event.b0_lp_reco_eta, event.b0_lp_reco_phi, eleM)
        lv_ele_sublead_reco.SetPtEtaPhiM(event.b0_lm_reco_pt, event.b0_lm_reco_eta, event.b0_lm_reco_phi, eleM)
        lv_ele_lead_gen.SetPtEtaPhiM(event.b0_lp_pt, event.b0_lp_eta, event.b0_lp_phi, eleM)
        lv_ele_sublead_gen.SetPtEtaPhiM(event.b0_lm_pt, event.b0_lm_eta, event.b0_lm_phi, eleM)
    else:
        lv_ele_lead_reco.SetPtEtaPhiM(event.b0_lm_reco_pt, event.b0_lm_reco_eta, event.b0_lm_reco_phi, eleM)
        lv_ele_sublead_reco.SetPtEtaPhiM(event.b0_lp_reco_pt, event.b0_lp_reco_eta, event.b0_lp_reco_phi, eleM)
        lv_ele_lead_gen.SetPtEtaPhiM(event.b0_lm_pt, event.b0_lm_eta, event.b0_lm_phi, eleM)
        lv_ele_sublead_gen.SetPtEtaPhiM(event.b0_lp_pt, event.b0_lp_eta, event.b0_lp_phi, eleM)

    if event.b0_k_pt > event.b0_pi_pt:
        lv_kaon_lead_gen.SetPtEtaPhiM(event.b0_k_pt, event.b0_k_eta, event.b0_k_phi, kaonM)
        lv_kaon_sublead_gen.SetPtEtaPhiM(event.b0_pi_pt, event.b0_pi_eta, event.b0_pi_phi, kaonM)
        lv_kaon_lead_reco.SetPtEtaPhiM(event.b0_k_reco_pt, event.b0_k_reco_eta, event.b0_k_reco_phi, kaonM)
        lv_kaon_sublead_reco.SetPtEtaPhiM(event.b0_pi_reco_pt, event.b0_pi_reco_eta, event.b0_pi_reco_phi, kaonM)
    else:
        lv_kaon_lead_gen.SetPtEtaPhiM(event.b0_pi_pt, event.b0_pi_eta, event.b0_pi_phi, kaonM)
        lv_kaon_sublead_gen.SetPtEtaPhiM(event.b0_k_pt, event.b0_k_eta, event.b0_k_phi, kaonM)
        lv_kaon_lead_reco.SetPtEtaPhiM(event.b0_pi_reco_pt, event.b0_pi_reco_eta, event.b0_pi_reco_phi, kaonM)
        lv_kaon_sublead_reco.SetPtEtaPhiM(event.b0_k_reco_pt, event.b0_k_reco_eta, event.b0_k_reco_phi, kaonM)

    lv_bs_reco = lv_ele_lead_reco + lv_ele_sublead_reco + lv_kaon_lead_reco + lv_kaon_sublead_reco

    #if (lv_ele_lead_reco + lv_ele_sublead_reco).M2() < 2: continue

    # Fill histograms

    h_q2_InvM_mc.Fill(pow((lv_ele_lead_reco + lv_ele_sublead_reco).M(),2))
    h_jpsi_InvM_mc.Fill((lv_ele_lead_reco + lv_ele_sublead_reco).M())
    h_phiee_InvM_mc.Fill((lv_kaon_lead_reco + lv_kaon_sublead_reco).M())
    h_bsee_InvM_mc.Fill((lv_ele_lead_reco + lv_ele_sublead_reco + lv_kaon_lead_reco + lv_kaon_sublead_reco).M())

    h_electron_pt_lead_mc.Fill(lv_ele_lead_reco.Pt())
    h_electron_eta_lead_mc.Fill(lv_ele_lead_reco.Eta())
    h_electron_phi_lead_mc.Fill(lv_ele_lead_reco.Phi())

    h_electron_pt_sublead_mc.Fill(lv_ele_sublead_reco.Pt())
    h_electron_eta_sublead_mc.Fill(lv_ele_sublead_reco.Eta())
    h_electron_phi_sublead_mc.Fill(lv_ele_sublead_reco.Phi())

    if event.b0_lp_reco_pt > event.b0_lm_reco_pt:
        h_electron_D0_lead_mc.Fill(event.b0_lp_reco_dxy)
        h_electron_Dz_lead_mc.Fill(event.b0_lp_reco_dz)
        h_electron_D0_sublead_mc.Fill(event.b0_lm_reco_dxy)
        h_electron_Dz_sublead_mc.Fill(event.b0_lm_reco_dz)

    else:
        h_electron_D0_lead_mc.Fill(event.b0_lm_reco_dxy)
        h_electron_Dz_lead_mc.Fill(event.b0_lm_reco_dz)
        h_electron_D0_sublead_mc.Fill(event.b0_lp_reco_dxy)
        h_electron_Dz_sublead_mc.Fill(event.b0_lp_reco_dz)

    if event.b0_k_reco_pt > event.b0_pi_reco_pt:
        h_kaon_ee_pt_lead_mc.Fill(event.b0_k_reco_pt)
        h_kaon_ee_eta_lead_mc.Fill(event.b0_k_reco_eta)
        h_kaon_ee_phi_lead_mc.Fill(event.b0_k_reco_phi)
        h_kaon_ee_Dz_lead_mc.Fill(event.b0_k_reco_dz)
        h_kaon_ee_D0_lead_mc.Fill(event.b0_k_reco_dxy)
        h_kaon_ee_D0Sig_lead_mc.Fill(event.b0_k_reco_dxy/event.b0_k_reco_dxyError)
        h_kaon_ee_normChi2_lead_mc.Fill(event.b0_k_reco_chi2/event.b0_k_reco_ndof if event.b0_k_reco_ndof != 0 else event.b0_k_reco_chi2*1.0e+6)

        h_kaon_ee_pt_sublead_mc.Fill(event.b0_pi_reco_pt)
        h_kaon_ee_eta_sublead_mc.Fill(event.b0_pi_reco_eta)
        h_kaon_ee_phi_sublead_mc.Fill(event.b0_pi_reco_phi)
        h_kaon_ee_Dz_sublead_mc.Fill(event.b0_pi_reco_dz)
        h_kaon_ee_D0_sublead_mc.Fill(event.b0_pi_reco_dxy)
        h_kaon_ee_D0Sig_sublead_mc.Fill(event.b0_pi_reco_dxy/event.b0_pi_reco_dxyError)
        h_kaon_ee_normChi2_sublead_mc.Fill(event.b0_pi_reco_chi2/event.b0_pi_reco_ndof if event.b0_pi_reco_ndof != 0 else event.b0_pi_reco_chi2*1.0e+6)

    else:
        h_kaon_ee_pt_lead_mc.Fill(event.b0_pi_reco_pt)
        h_kaon_ee_eta_lead_mc.Fill(event.b0_pi_reco_eta)
        h_kaon_ee_phi_lead_mc.Fill(event.b0_pi_reco_phi)
        h_kaon_ee_Dz_lead_mc.Fill(event.b0_pi_reco_dz)
        h_kaon_ee_D0_lead_mc.Fill(event.b0_pi_reco_dxy)
        h_kaon_ee_D0Sig_lead_mc.Fill(event.b0_pi_reco_dxy/event.b0_pi_reco_dxyError)
        h_kaon_ee_normChi2_lead_mc.Fill(event.b0_pi_reco_chi2/event.b0_pi_reco_ndof if event.b0_pi_reco_ndof != 0 else event.b0_pi_reco_chi2*1.0e+6)

        h_kaon_ee_pt_sublead_mc.Fill(event.b0_k_reco_pt)
        h_kaon_ee_eta_sublead_mc.Fill(event.b0_k_reco_eta)
        h_kaon_ee_phi_sublead_mc.Fill(event.b0_k_reco_phi)
        h_kaon_ee_Dz_sublead_mc.Fill(event.b0_k_reco_dz)
        h_kaon_ee_D0_sublead_mc.Fill(event.b0_k_reco_dxy)
        h_kaon_ee_D0Sig_sublead_mc.Fill(event.b0_k_reco_dxy/event.b0_k_reco_dxyError)
        h_kaon_ee_normChi2_sublead_mc.Fill(event.b0_k_reco_chi2/event.b0_k_reco_ndof if event.b0_k_reco_ndof != 0 else event.b0_k_reco_chi2*1.0e+6)


    ctxy = ((event.b0_reco_x - event.b0_reco_vtx)*lv_bs_reco.Px() + (event.b0_reco_y - event.b0_reco_vty)*lv_bs_reco.Py())/(lv_bs_reco.Pt()**2)*bsM


    h_bs_pt_ee_mc.Fill(lv_bs_reco.Pt())
    h_svProb_ee_mc.Fill(ROOT.TMath.Prob(event.b0_reco_chi2, int(event.b0_reco_ndof)))  
    h_svCosAngle_ee_mc.Fill(event.b0_reco_cosAngle)
    h_svCtxy_ee_mc.Fill(ctxy)
    h_svLxy_ee_mc.Fill(event.b0_reco_Lxy)
    h_svLxySig_ee_mc.Fill(event.b0_reco_Lxy/event.b0_reco_LxyError)    
    
    h_electron_dR_mc.Fill(lv_ele_lead_reco.DeltaR(lv_ele_sublead_reco))
    h_kaon_ee_dR_mc.Fill(lv_kaon_lead_reco.DeltaR(lv_kaon_sublead_reco))
    h_jpsiPhiOpen_ee_mc.Fill((lv_ele_lead_reco + lv_ele_sublead_reco).DeltaR(lv_kaon_lead_reco + lv_kaon_sublead_reco))

    nTot = nTot + 1

print("Total number of events: " + str(nTot))

file_out.cd()

file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0,2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0,2)) + ' minutes')
