#! /usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-s", "--signal", dest="signal", default="DoubleMuonNtu_Run2016B.list", help="List of MC siognal files")
parser.add_argument("-b", "--background", dest="background", default="DoubleMuonNtu_Run2016B.list", help="List of background input ggNtuplizer files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-m", "--maxevents", dest="maxevents", type=int, default=50000, help="Maximum number events to loop over")
parser.add_argument("-t", "--ttree", dest="ttree", default="ggNtuplizer/EventTree", help="TTree Name")
parser.add_argument("-v", "--trigger", dest="trigger", default=False, action='store_true', help="Boolean for choosing trigger or non-trigger events")
args = parser.parse_args()

import numpy as np
import itertools
import ROOT
import os
from collections import OrderedDict
from array import array
import math

ROOT.gROOT.ProcessLine(
"struct bsFeatures_t {\
  Float_t       elePtLead;\
  Float_t       eleEtaLead;\
  Float_t       elePhiLead;\
  Float_t       eleD0Lead;\
  Float_t       eleDzLead;\
  Float_t       elePtSublead;\
  Float_t       eleEtaSublead;\
  Float_t       elePhiSublead;\
  Float_t       eleD0Sublead;\
  Float_t       eleDzSublead;\
  Float_t       kaonPtLead;\
  Float_t       kaonEtaLead;\
  Float_t       kaonPhiLead;\
  Float_t       kaonD0Lead;\
  Float_t       kaonDzLead;\
  Float_t       kaonNormChi2Lead;\
  Float_t       kaonPtSublead;\
  Float_t       kaonEtaSublead;\
  Float_t       kaonPhiSublead;\
  Float_t       kaonD0Sublead;\
  Float_t       kaonDzSublead;\
  Float_t       kaonNormChi2Sublead;\
  Float_t       eledR;\
  Float_t       kaondR;\
  Float_t       jpsiPhidR;\
  Float_t       svCtxy;\
  Float_t       svProb;\
  Float_t       svCosine;\
  Float_t       svLxy;\
  Float_t       svLxySig;\
  Float_t       jpsiPt;\
  Float_t       phiPt;\
  Float_t       bsPt;\
  Float_t       jpsiMass;\
  Float_t       phiMass;\
  Float_t       bsMass;\
};" );

if os.path.isfile('~/.rootlogon.C'): ROOT.gROOT.Macro(os.path.expanduser('~/.rootlogon.C'))
ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(000000)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.UseCurrentStyle()

sw = ROOT.TStopwatch()
sw.Start()

eleM = 0.0005109989461
muonM = 0.1056583745
kaonM = 0.493677
bsM = 5.3663
bssideminuslow = 5.1
bssideminusup = 5.3
bssidepluslow = 5.45
bssideplusup = 5.6
jpsisideminuslow = 2.4
jpsisideminusup = 2.8
jpsisidepluslow = 3.3
jpsisideplusup = 3.8
phisideminuslow = 1.0
phisideminusup = 1.01
phisidepluslow = 1.04
phisideplusup = 1.05
bslow = 4.8
bsup = 6.0

lv_ele_lead_reco = ROOT.TLorentzVector()
lv_ele_sublead_reco = ROOT.TLorentzVector()
lv_kaon_lead_reco = ROOT.TLorentzVector()
lv_kaon_sublead_reco = ROOT.TLorentzVector()
lv_bs_reco = ROOT.TLorentzVector()

# Output file and any histograms we want
file_out = ROOT.TFile(args.outputfile, 'recreate')

branches = [
  'elePtLead',
  'eleEtaLead',
  'elePhiLead',
  'eleD0Lead',
  'eleDzLead',
  'elePtSublead',
  'eleEtaSublead',
  'elePhiSublead',
  'eleD0Sublead',
  'eleDzSublead',
  'kaonPtLead',
  'kaonEtaLead',
  'kaonPhiLead',
  'kaonD0Lead',
  'kaonDzLead',
  'kaonNormChi2Lead',
  'kaonPtSublead',
  'kaonEtaSublead',
  'kaonPhiSublead',
  'kaonD0Sublead',
  'kaonDzSublead',
  'kaonNormChi2Sublead',
  'eledR',
  'kaondR',
  'jpsiPhidR',
  'svCtxy',
  'svProb',
  'svCosine',
  'svLxy',
  'svLxySig',
  'jpsiPt',
  'phiPt',
  'bsPt',
  'jpsiMass',
  'phiMass',
  'bsMass'
]

#branch_descriptor = ':'.join([br + '/F' for br in branches])

print args.maxevents, args.trigger

numBGEvents = 0
numSGEvents = 0
foundBGEvent = False
rand = ROOT.TRandom3()
rand.SetSeed(0)

MCtchain = ROOT.TChain('tree')
MCtchain.Add(args.signal)
print 'Total number of signal events: ' + str(MCtchain.GetEntries())

bsFeatures = ROOT.bsFeatures_t()
stree = ROOT.TTree('signal', 'signal tree')
for var in branches:
  stree.Branch(var, ROOT.AddressOf(bsFeatures, var), var+'/F')


for ievent,event in enumerate(MCtchain):
  #if ievent > args.maxevents and args.maxevents != -1: break
  #if ievent % 1000 == 0: print 'Processing entry ' + str(ievent)
  if ievent % 10000 == 0: print 'Processing entry ' + str(ievent) +' , number of selected events: ' + str(numSGEvents)
  if numSGEvents >= args.maxevents: break
  
  if event.tag_mu_reco_dxy == -99: continue
  if abs(event.b0_lp_eta) > 2.5 or abs(event.b0_lm_eta) > 2.5 or abs(event.b0_k_eta) > 2.5 or abs(event.b0_pi_eta) > 2.5: continue
  if abs(event.b0_lp_reco_pt) < 0.4 or abs(event.b0_lm_reco_pt) < 0.4 or abs(event.b0_k_reco_pt) < 0.4 or abs(event.b0_pi_reco_pt) < 0.4: continue
  if event.foundB0 < 0.7: continue
  #if event.tag_mu_reco_pt < 9: continue
  #if abs(event.tag_mu_reco_dxy/event.tag_mu_reco_dxyError) < 6: continue

  if event.b0_lp_pt > event.b0_lm_pt:
      lv_ele_lead_reco.SetPtEtaPhiM(event.b0_lp_reco_pt, event.b0_lp_reco_eta, event.b0_lp_reco_phi, eleM)
      lv_ele_sublead_reco.SetPtEtaPhiM(event.b0_lm_reco_pt, event.b0_lm_reco_eta, event.b0_lm_reco_phi, eleM)
  else:
      lv_ele_lead_reco.SetPtEtaPhiM(event.b0_lm_reco_pt, event.b0_lm_reco_eta, event.b0_lm_reco_phi, eleM)
      lv_ele_sublead_reco.SetPtEtaPhiM(event.b0_lp_reco_pt, event.b0_lp_reco_eta, event.b0_lp_reco_phi, eleM)

  lv_kaon_lead_reco.SetPtEtaPhiM(event.b0_k_reco_pt, event.b0_k_reco_eta, event.b0_k_reco_phi, kaonM)
  lv_kaon_sublead_reco.SetPtEtaPhiM(event.b0_pi_reco_pt, event.b0_pi_reco_eta, event.b0_pi_reco_phi, kaonM)

  lv_bs_reco = lv_ele_lead_reco + lv_ele_sublead_reco + lv_kaon_lead_reco + lv_kaon_sublead_reco

  bsFeatures.elePtLead = lv_ele_lead_reco.Pt()/lv_bs_reco.M()
  bsFeatures.eleEtaLead = lv_ele_lead_reco.Eta()
  bsFeatures.elePhiLead = lv_ele_lead_reco.Phi()
  
  bsFeatures.elePtSublead = lv_ele_sublead_reco.Pt()/lv_bs_reco.M()
  bsFeatures.eleEtaSublead = lv_ele_sublead_reco.Eta()
  bsFeatures.elePhiSublead = lv_ele_sublead_reco.Phi()

  if event.b0_lp_reco_pt > event.b0_lm_reco_pt:
    bsFeatures.eleD0Lead = event.b0_lp_reco_dxy
    bsFeatures.eleDzLead = event.b0_lp_reco_dz
    bsFeatures.eleD0Sublead = event.b0_lm_reco_dxy
    bsFeatures.eleDzSublead = event.b0_lm_reco_dz

  else:
    bsFeatures.eleD0Lead = event.b0_lm_reco_dxy
    bsFeatures.eleDzLead = event.b0_lm_reco_dz
    bsFeatures.eleD0Sublead = event.b0_lp_reco_dxy
    bsFeatures.eleDzSublead = event.b0_lp_reco_dz

  if event.b0_k_reco_pt > event.b0_pi_reco_pt:
    bsFeatures.kaonPtLead = event.b0_k_reco_pt/lv_bs_reco.M()
    bsFeatures.kaonEtaLead = event.b0_k_reco_eta
    bsFeatures.kaonPhiLead = event.b0_k_reco_phi
    bsFeatures.kaonD0Lead = event.b0_k_reco_dxy
    bsFeatures.kaonDzLead = event.b0_k_reco_dz
    bsFeatures.kaonNormChi2Lead = event.b0_k_reco_chi2/event.b0_k_reco_ndof if event.b0_k_reco_ndof != 0 else event.b0_k_reco_chi2*1.0e+6

    bsFeatures.kaonPtSublead = event.b0_pi_reco_pt/lv_bs_reco.M()
    bsFeatures.kaonEtaSublead = event.b0_pi_reco_eta
    bsFeatures.kaonPhiSublead = event.b0_pi_reco_phi
    bsFeatures.kaonD0Sublead = event.b0_pi_reco_dxy
    bsFeatures.kaonDzSublead = event.b0_pi_reco_dz
    bsFeatures.kaonNormChi2Sublead = event.b0_pi_reco_chi2/event.b0_pi_reco_ndof if event.b0_pi_reco_ndof != 0 else event.b0_pi_reco_chi2*1.0e+6

  else:
    bsFeatures.kaonPtLead = event.b0_pi_reco_pt/lv_bs_reco.M()
    bsFeatures.kaonEtaLead = event.b0_pi_reco_eta
    bsFeatures.kaonPhiLead = event.b0_pi_reco_phi
    bsFeatures.kaonD0Lead = event.b0_pi_reco_dxy
    bsFeatures.kaonDzLead = event.b0_pi_reco_dz
    bsFeatures.kaonNormChi2Lead = event.b0_pi_reco_chi2/event.b0_pi_reco_ndof if event.b0_pi_reco_ndof != 0 else event.b0_pi_reco_chi2*1.0e+6

    bsFeatures.kaonPtSublead = event.b0_k_reco_pt/lv_bs_reco.M()
    bsFeatures.kaonEtaSublead = event.b0_k_reco_eta
    bsFeatures.kaonPhiSublead = event.b0_k_reco_phi
    bsFeatures.kaonD0Sublead = event.b0_k_reco_dxy
    bsFeatures.kaonDzSublead = event.b0_k_reco_dz
    bsFeatures.kaonNormChi2Sublead = event.b0_k_reco_chi2/event.b0_k_reco_ndof if event.b0_k_reco_ndof != 0 else event.b0_k_reco_chi2*1.0e+6

  bsFeatures.eledR = lv_ele_lead_reco.DeltaR(lv_ele_sublead_reco)
  bsFeatures.kaondR = lv_kaon_lead_reco.DeltaR(lv_kaon_sublead_reco)
  bsFeatures.jpsiPhidR = (lv_kaon_lead_reco + lv_kaon_sublead_reco).DeltaR(lv_ele_lead_reco + lv_ele_sublead_reco)

  bsFeatures.jpsiPt = (lv_ele_lead_reco + lv_ele_sublead_reco).Pt()/lv_bs_reco.M()
  bsFeatures.phiPt = (lv_kaon_lead_reco + lv_kaon_sublead_reco).Pt()/lv_bs_reco.M()
  bsFeatures.bsPt = lv_bs_reco.Pt()/lv_bs_reco.M()

  ctxy = ((event.b0_reco_x - event.b0_reco_vtx)*lv_bs_reco.Px() + (event.b0_reco_y - event.b0_reco_vty)*lv_bs_reco.Py())/(lv_bs_reco.Pt()**2)*bsM

  bsFeatures.svCtxy = ctxy
  bsFeatures.svProb = ROOT.TMath.Prob(event.b0_reco_chi2, int(event.b0_reco_ndof))
  bsFeatures.svCosine = event.b0_reco_cosAngle
  bsFeatures.svLxy = event.b0_reco_Lxy
  bsFeatures.svLxySig = event.b0_reco_Lxy/event.b0_reco_LxyError

  bsFeatures.jpsiMass = (lv_ele_lead_reco + lv_ele_sublead_reco).M()
  bsFeatures.phiMass = (lv_kaon_lead_reco + lv_kaon_sublead_reco).M()
  bsFeatures.bsMass = lv_bs_reco.M()

  stree.Fill()
  numSGEvents = numSGEvents + 1



file_out.cd()

file_out.Write()
file_out.Close()

sw.Stop()
print 'Real time: ' + str(round(sw.RealTime() / 60.0,2)) + ' minutes'
print 'CPU time:  ' + str(round(sw.CpuTime() / 60.0,2)) + ' minutes'
