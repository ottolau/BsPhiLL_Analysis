import os
import sys
import ROOT
from itertools import combinations
import numpy as np
from array import array


ROOT.gROOT.SetBatch(ROOT.kTRUE);
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)

varList = ["pt",
           "eta",
           "phi",
           "dR",
           "D0",
           "Dz",
           "SIP",
           "D0Error",
           "DzError",
           "D0Sig",
           "DzSig",
           "normChi2",
           "svCtxy",
           "svChi2",
           "svLxy",
           "svProb",
           "svCosAngle",
           "electron_dR",
           "kaon_ee_dR",
           "jpsiPhiOpen",
           "jpsi_InvM",
           "phiee_InvM",
           "bsee_InvM",
           "q2"
            ]

varUnitMap = {"pt": "p_{T} [GeV]",
              "eta": "#eta",
              "phi": "#phi",
              "dR": "#Delta R",
              "D0": "d_{xy} [cm]",
              "Dz": "d_{z} [cm]",
              "SIP": "dB/#sigma_{dB}",
              "D0Error": "#sigma_{d_{xy}} [cm]",
              "DzError": "#sigma_{d_{z}} [cm]",
              "D0Sig": "d_{xy}/#sigma_{d_{xy}}",
              "DzSig": "d_{z}/#sigma_{d_{z}}",
              "normChi2": "#chi^{2}_{track}/d.o.f.",
              "svCtxy": "ct_{xy} [cm]",
              "svLxy": "L_{xy} / #sigma_{L_{xy}}",
              "svChi2": "#chi^{2}_{SV}",
              "svProb": "P(#chi^{2}_{SV})",
              "svCosAngle": "cos #alpha_{2D}",
              "electron_dR": "#Delta R(e^{+}, e^{-})",
              "kaon_ee_dR": "#Delta R(K^{+}, K^{-})",
              "jpsiPhiOpen": "#Delta R(e^{+}e^{-}, #phi)",
              "jpsi_InvM": "m(e^{+}e^{-}) [GeV]",
              "phiee_InvM": "m(K^{+}K^{-}) [GeV]",
              "bsee_InvM": "m(K^{+}K^{-}e^{+}e^{-}) [GeV]",
              "q2": "q^{2} [GeV^{2}]"
                }


def make_plots(filename):
    f = ROOT.TFile(filename)
    dir_list = ROOT.gDirectory.GetListOfKeys()
    outputfile = "BsPhiJpsiEE_KinematicsPlots_mc"
    c_list = []

    for key in dir_list:
        if key.GetClassName() == "TH1D":
            histo = key.ReadObj()
            histo_name = histo.GetName()
            canvas_name = "c_" + histo_name
            c = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
            unit = ""
            for var in varList:
                if var in histo_name: unit = varUnitMap[var]
           
            histo.SetTitle(histo_name)
            histo.GetYaxis().SetTitle("N")
            histo.GetXaxis().SetTitle(unit)
            histo.SetLineColor(2)
            histo.DrawCopy()             
            c_list.append(c)

    for i,c in enumerate(c_list):
        if i == 0:
            c.Print("Figures/{}.pdf(".format(outputfile),"pdf")
        elif i == len(c_list)-1:
            c.Print("Figures/{}.pdf)".format(outputfile),"pdf")
        else:
            c.Print("Figures/{}.pdf".format(outputfile),"pdf")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Produce histograms')
    parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="", help="ROOT file contains histograms")
    args = parser.parse_args()

    make_plots(args.inputfiles)

