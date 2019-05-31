#! /usr/bin/env python

import argparse

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="DoubleMuonNtu_Run2016B.list", help="List of input ggNtuplizer files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-m", "--maxevents", dest="maxevents", type=int, default=-1, help="Maximum number events to loop over")
parser.add_argument("-t", "--ttree", dest="ttree", default="cutTree", help="TTree Name")
parser.add_argument("-c", "--mvacut", dest="mvacut", type=float, default=0.5, help="MVA cut")
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

h_q2_InvM = ROOT.TH1D("h_q2_InvM", "", 100, 0.0, 25.0);
h_jpsi_InvM = ROOT.TH1D("h_jpsi_InvM", "", 100, 0.0, 5.0);
h_phiee_InvM = ROOT.TH1D("h_phiee_InvM", "", 50, 0.98, 1.06);
h_bsee_InvM = ROOT.TH1D("h_bsee_InvM", "", 30, 4.5, 6.0);
h_bsee_InvM_lowerPhiSB = ROOT.TH1D("h_bsee_InvM_lowerPhiSB", "", 30, 4.5, 6.0);
h_bsee_InvM_upperPhiSB = ROOT.TH1D("h_bsee_InvM_upperPhiSB", "", 30, 4.5, 6.0);
#h_bsee_InvM = ROOT.TH1D("h_bsee_InvM", "", 20, 4.8, 6.0);
#h_elePt_lead = ROOT.TH1D("h_elePt_lead", "", 40, 0.0, 20.0);
#h_elePt_sublead = ROOT.TH1D("h_elePt_sublead", "", 40, 0.0, 20.0);
#h_kaonPt_lead = ROOT.TH1D("h_kaonPt_lead", "", 40, 0.0, 20.0);
#h_kaonPt_sublead = ROOT.TH1D("h_kaonPt_sublead", "", 40, 0.0, 20.0);
#h_electron_dR = ROOT.TH1D("h_electron_dR", "", 120, 0.0, 4.0);
#h_kaon_ee_dR = ROOT.TH1D("h_kaon_ee_dR", "", 120, 0.0, 4.0);
#h_jpsiPhiOpen_ee = ROOT.TH1D("h_jpsiPhiOpen_ee", "", 120, 0.0, 4.0);


eleM = 0.0005109989461
kaonM = 0.493677
pionM = 0.13957018
bsM = 5.3663
jpsilow = 2.9
jpsiup = 3.25
#jpsilow = 2.5
#jpsiup = 2.8
philow = 1.015
phiup = 1.03
#philow = 1.04
#phiup = 1.06
#philow = 0.98
#phiup = 1.01
blow = 5.0
bup = 5.5

bslow = 5.2
bsup = 5.6

temp_m_ee = 0.

#Loop over all the events in the input ntuple
for ievent,event in enumerate(tchain):
    if ievent > args.maxevents and args.maxevents != -1: break
    if ievent % 100000 == 0: print('Processing entry ' + str(ievent))

    #if event.kaonPtLead < 2.0: continue
    if event.bsPt < 4.0: continue
    if event.LxySig < 2.0: continue
    if event.SvProb < 0.1: continue
    if event.SvCosine < 0.999: continue
    #if event.mvaUnBWP_sublead < 4.5: continue
    if event.elePtLead < 2.0 or event.elePtSublead < 2.0: continue
    #if event.elePtLead > 2.0 and event.elePtSublead > 2.0: continue
    if event.mvaUnBWP_lead < 3.05 or event.mvaUnBWP_sublead < 3.05: continue
    #if event.mvaUnBWP_lead < 5.0 or event.mvaUnBWP_sublead < 5.0: continue
    #if event.mvaUnBWP_lead < 6.0 or event.mvaUnBWP_sublead < 0.19: continue

    if event.m_ee != temp_m_ee:
        h_q2_InvM.Fill(pow(event.m_ee, 2))
        h_jpsi_InvM.Fill(event.m_ee)
        temp_m_ee = event.m_ee
    h_phiee_InvM.Fill(event.m_KK)
    #if jpsilow < event.m_ee < jpsiup and philow < event.m_KK < phiup:
    #if philow < event.m_KK < phiup:
    #if jpsilow < event.m_ee < jpsiup:
    #if jpsilow < event.m_ee < jpsiup and philow < event.m_KK < phiup and event.eledR < 1.5 and event.kaondR < 0.5 and event.jpsiPhidR < 1.5:
    #if jpsilow < event.m_ee < jpsiup and philow < event.m_KK < phiup and event.eledR > 0.3:
    if jpsilow < event.m_ee < jpsiup:
        if philow < event.m_KK < phiup:
            h_bsee_InvM.Fill(event.m_KKee)
        if 0.98 < event.m_KK < 1.01:
            h_bsee_InvM_lowerPhiSB.Fill(event.m_KKee)
        if 1.04 < event.m_KK < 1.06:
            h_bsee_InvM_upperPhiSB.Fill(event.m_KKee)
        #if blow < event.m_KKee < bup:
        #if bslow < event.m_KKee < bsup:
            #h_elePt_lead.Fill(event.elePtLead)
            #h_elePt_sublead.Fill(event.elePtSublead)
            #h_kaonPt_lead.Fill(event.kaonPtLead)
            #h_kaonPt_sublead.Fill(event.kaonPtSublead)
    #h_electron_dR.Fill(event.eledR)
    #h_kaon_ee_dR.Fill(event.kaondR)
    #h_jpsiPhiOpen_ee.Fill(event.jpsiPhidR)


file_out.cd()

file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0,2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0,2)) + ' minutes')
