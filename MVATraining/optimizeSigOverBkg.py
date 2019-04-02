#!/usr/bin/env python
import os,sys
from time import sleep
import math
import numpy as np
import ROOT
from ROOT import RooFit
#-------------------------------------------------------------
eleM = 0.0005109989461
kaonM = 0.493677
pionM = 0.13957018
bsM = 5.3663
jpsilow = 2.9
jpsiup = 3.2
philow = 1.01
phiup = 1.03
#-------------------------------------------------------------

def fitHisto(tchain, outputfile, mvaCut, index, mvaList):
    hist = ROOT.TH1D("hist", "", 25, 4.5, 6.0);
    #Loop over all the events in the input ntuple
    print("Filling and fitting mva cut at %f..."%(mvaCut))
    for ievent,event in enumerate(tchain):
        #if ievent % 100000 == 0: print('Processing entry ' + str(ievent))

        if event.mvaValue < mvaCut: continue
        if jpsilow < event.m_ee < jpsiup and philow < event.m_KK < phiup:
            hist.Fill(event.m_KKee)

    #hist.Rebin(2)
    hist.SetStats(0)
    hdata = ROOT.RooDataHist("hdata", "hdata", l, RooFit.Import(hist))

    results = model.fitTo(hdata, RooFit.Save(), RooFit.Extended(True))
    #results.Print()

    #bkgRange = ROOT.RooRealVar("bkgRange","bkgRange",5.25,5.45) ;
    #bkgRangeArgSet = ROOT.RooArgSet(bkgRange)
    #fracBkgRange = bkg.createIntegral(bkgRangeArgSet,"window") ;

    x.setRange("window",5.1,5.5) 
    fracBkgRange = bkg.createIntegral(obs,obs,"window") 
    nbkgWindow = nbkg.getVal() * fracBkgRange.getVal() 
    nbkgWindowError = fracBkgRange.getPropagatedError(results) * nbkg.getVal()
    print("MAV cut: %f, Number of signals: %f, Number of background: %f, S/sqrt(S+B): %f"%(mvaCut, nsig.getVal(), nbkgWindow, nsig.getVal()/math.sqrt(nsig.getVal() + nbkgWindow)))

    # Plot results of fit on a different frame
    c2 = ROOT.TCanvas('fig_binnedFit', 'fit', 700, 500)
    c2.SetGrid()
    c2.cd()
    ROOT.gPad.SetLeftMargin(0.10)
    ROOT.gPad.SetRightMargin(0.05)

    xframe = wspace.var('x').frame(RooFit.Title("MVA Cut: %f"%(mvaCut)))
    hdata.plotOn(xframe, RooFit.Name("data"))
    model.plotOn(xframe,RooFit.Name("global"),RooFit.LineColor(2),RooFit.MoveToBack()) # this will show fit overlay on canvas
    model.plotOn(xframe,RooFit.Name("bkg"),RooFit.Components("bkg"),RooFit.LineStyle(ROOT.kDashed),RooFit.LineColor(ROOT.kMagenta),RooFit.MoveToBack()) ;
    model.plotOn(xframe,RooFit.Name("sig"),RooFit.Components("sig"),RooFit.DrawOption("FL"),RooFit.FillColor(9),RooFit.FillStyle(3004),RooFit.LineStyle(6),RooFit.LineColor(9)) ;
    model.plotOn(xframe,RooFit.VisualizeError(results), RooFit.FillColor(ROOT.kOrange), RooFit.MoveToBack()) # this will show fit overlay on canvas

    xframe.GetXaxis().SetTitleFont(42)
    xframe.GetXaxis().SetLabelFont(42)
    xframe.GetYaxis().SetTitleFont(42)
    xframe.GetYaxis().SetLabelFont(42)
    xframe.GetXaxis().SetTitle("m(K^{+}K^{-}e^{+}e^{-}) [GeV]")
    xframe.SetStats(0)
    xframe.SetMinimum(0)
    xframe.Draw()

    legend = ROOT.TLegend(0.65,0.65,0.9,0.85);
    legend.SetTextFont(42);
    legend.SetTextSize(0.04);
    legend.AddEntry(xframe.findObject("data"),"Data","lpe");
    legend.AddEntry(xframe.findObject("bkg"),"Background fit","l");
    legend.AddEntry(xframe.findObject("sig"),"Signal fit","l");
    legend.AddEntry(xframe.findObject("global"),"Global Fit","l");
    legend.Draw();

    #c2.SaveAs('%s_mvaCut_%s.png'%(outputfile, str(mvaCut).replace('.','_')))
    if index == 0:
        c2.Print("{}.pdf(".format(outputfile),"pdf")
    elif index == len(mvaList)-1:
        c2.Print("{}.pdf)".format(outputfile),"pdf")
    else:
        c2.Print("{}.pdf".format(outputfile),"pdf")

    return nsig.getVal(), nsig.getError(), nbkgWindow, nbkgWindowError, nsig.getVal()/math.sqrt(nsig.getVal() + nbkgWindow)

#-------------------------------------------------------------

def plotSNR(outputfile, cut, sig, sigError, SNR):
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib import rc
    #.Allow for using TeX mode in matplotlib Figures
    rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
    rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]

    ratio=5.0/7.0
    fig_width_pt = 3*246.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = ratio if ratio != 0.0 else (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]

    params = {'text.usetex' : True,
            'axes.labelsize': 24,
            'font.size': 24,
            'legend.fontsize': 20,
            'xtick.labelsize': 24,
            'ytick.labelsize': 24,
            'font.family' : 'lmodern',
            'text.latex.unicode': True,
            'axes.grid' : False,
            'text.usetex': True,
            'figure.figsize': fig_size}
    plt.rcParams.update(params)

    fig, ax1 = plt.subplots()
    plt.grid(linestyle='--')
    ax1.plot(cut, SNR, 'b-', label=r'$S/\sqrt{S+B}$')
    ax1.set_xlabel(r'MVA Cut')
    ax1.set_ylabel(r'$S/\sqrt{S+B}$')
    #ax1.set_ylim(ymin=0)
    #for t1 in ax1.get_yticklabels():
        #t1.set_color('b')

    lower_bound = [s-serror for (s, serror) in zip(sig, sigError)]
    upper_bound = [s+serror for (s, serror) in zip(sig, sigError)]

    ax2 = ax1.twinx()
    ax2.plot(cut, sig, 'r--', label=r'Number of signals')
    ax2.fill_between(cut, lower_bound, upper_bound, facecolor='yellow', alpha=0.5)
    ax2.set_ylabel(r'Number of signals')
    ax2.set_ylim(ymin=0)
    #for t1 in ax2.get_yticklabels():
        #t1.set_color('r')
    #ax2.legend(loc='upper left')
    #lns = lns1+lns2
    #labs = [l.get_label() for l in lns]
    #ax2.legend(lns, labs, loc='upper left')
    #fig.legend(loc=2, bbox_to_anchor=(0,1), bbox_transform=ax2.transAxes)
    fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=ax2.transAxes)
    fig.savefig('%s_SNRPlot.pdf'%(outputfile), bbox_inches='tight')


#msgservice = ROOT.RooMsgService.instance()
#msgservice.setGlobalKillBelow(RooFit.FATAL)
wspace = ROOT.RooWorkspace('myWorkSpace')
ROOT.gStyle.SetOptFit(0000);
ROOT.gROOT.SetBatch(True);
ROOT.gROOT.SetStyle("Plain");
ROOT.gStyle.SetGridStyle(3);
ROOT.gStyle.SetOptStat(000000);


inputfile = '/eos/uscms/store/user/klau/BsPhiLL_output/BsPhiEE_MVAEvaluation/BsPhiJpsiEE_MVATree_noPhiM_BOBest_samesign.root'
outputfile = 'MVACutOptimization/BsPhiJpsiEE_noPhiM_BOBest_samesign/MVAOpt_BsPhiJpsiEE_noPhiM_BOBest_samesign'
tchain = ROOT.TChain('mvaTree')
tchain.Add(inputfile)

xmin = 4.5
xmax = 6.0 
wspace.factory('x[5.0,%f,%f]' % (xmin, xmax))

wspace.factory('c0[1.0, -1.0, 1.0]')
wspace.factory('c1[-0.1, -1.0, 1.0]')
wspace.factory('c2[-0.1, -1.0, 1.0]')
wspace.factory('mean[5.3240e+00, 5.3240e+00, 5.3240e+00]')
wspace.factory('width[1.000e-01, 1.000e-01, 1.000e-01]')
wspace.factory('sigma[7.1858e-02, 7.1858e-02, 7.1858e-02]')
wspace.factory('nsig[100.0, 0.0, 100000.0]')
wspace.factory('nbkg[500.0, 0.0, 100000.0]')
wspace.factory('alpha[-1.0, -100.0, 0.0]')

x = wspace.var('x')
c0 = wspace.var('c0')
c1 = wspace.var('c1')
c2 = wspace.var('c2')
mean = wspace.var('mean')
width = wspace.var('width')
sigma = wspace.var('sigma')
nsig = wspace.var('nsig')
nbkg = wspace.var('nbkg')
alpha = wspace.var('alpha')

parameters = ['c0', 'c1', 'c2', 'mean', 'width', 'sigma', 'nsig', 'nbkg']

# NUMBER OF PARAMETERS
P = len(parameters)

wspace.factory('Voigtian::sig(x,mean,width,sigma)')
#wspace.factory('Chebychev::bkg(x,{c0,c1,c2})')
wspace.factory('Exponential::bkg(x,alpha)')
wspace.factory('SUM::model(nsig*sig,nbkg*bkg)')
         
model = wspace.pdf('model')
bkg = wspace.pdf('bkg')
sig = wspace.pdf('sig')

wspace.defineSet('obs', 'x')

# make the set obs known to Python
obs  = wspace.set('obs')
l = ROOT.RooArgList(obs)

SigList = []
SErrorList = []
BkgList = []
BErrorList = []
SNRList = []

mvaValue = np.linspace(0.5, 0.999, 100)
for i, mva in enumerate(mvaValue):
    numSig, sigError, numBkg, bkgError, SNR = fitHisto(tchain, outputfile, mva, i, mvaValue)
    SigList.append(numSig)
    SErrorList.append(sigError)
    BkgList.append(numBkg)
    BErrorList.append(bkgError)
    SNRList.append(SNR)

plotSNR(outputfile, mvaValue, SigList, SErrorList, SNRList)

