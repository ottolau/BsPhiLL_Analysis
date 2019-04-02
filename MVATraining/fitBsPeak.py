#!/usr/bin/env python
# Load the operating system and system modules into memory
import os,sys

# Load the sleep function from the time module
from time import sleep

# Import all functions from the math module if you are sure there
# will be no name collisions
import math

# Load everything that is in PyROOT. This is fine if, again, you are
# sure there are no name collisions between PyROOT and the names
# in other modules, otherwise use the syntax:
# example:
#  from ROOT import TCanvas, TFile, kBlue
#
import ROOT
from ROOT import RooFit
#-------------------------------------------------------------
# The main function can be called whatever you like. We follow
# C++ and call it main..dull, but clear!
def main():

    #msgservice = ROOT.RooMsgService.instance()
    #msgservice.setGlobalKillBelow(RooFit.FATAL)
    wspace = ROOT.RooWorkspace('myWorkSpace')
    ROOT.gStyle.SetOptFit(0000);
    ROOT.gROOT.SetBatch(True);
    ROOT.gROOT.SetStyle("Plain");
    ROOT.gStyle.SetGridStyle(3);
    ROOT.gStyle.SetOptStat(000000);

    xmin = 4.5 
    xmax = 6.0 
    
    wspace.factory('x[5.0,%f,%f]' % (xmin, xmax))

    #M = 15 
    #wspace.var('x').setBins(M)

    wspace.factory('c0[1.0, -1.0, 1.0]')
    wspace.factory('c1[-0.1, -1.0, 1.0]')
    wspace.factory('c2[-0.1, -1.0, 1.0]')
    wspace.factory('mean[5.36635e+00, 5.2, 5.6]')
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
    #wspace.factory("fsig[0.5,0.,1.]") 
    #model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
    #wspace.factory('AddPdf::model({sig,bkg},fsig)') 
    #wspace.factory('AddPdf::model({sig,bkg}, {nsig,nbkg})')
            
    model = wspace.pdf('model')
    #model = wspace.pdf('Voigtian')
    bkg = wspace.pdf('bkg')
    sig = wspace.pdf('sig')

    # define the set obs = (x)
    wspace.defineSet('obs', 'x')

    # make the set obs known to Python
    obs  = wspace.set('obs')
    l = ROOT.RooArgList(obs)

    f = ROOT.TFile("test.root","READ")
    hist = f.Get("h_bsee_InvM")
    hist.Rebin(2)
    hist.SetStats(0)
    hdata = ROOT.RooDataHist("hdata", "hdata", l, RooFit.Import(hist))

    #hdata = RooDataHist('hdata', 'binned data', obs)
    #hdata.add(data)  # add the data to the RooDataHist and bin them
    print("="*40)
    hdata.Print('verbose')
    print("="*40)

    # Do a multinomial fit to the binned data by
    # turning off extended likelihood mode. If you
    # want a multi-Poisson fit, change False to True.
    # (If interested, ask what all this means!)
    results = model.fitTo(hdata, RooFit.Save(), RooFit.Extended(True))
    #results = model.fitTo(hdata, RooFit.Minimizer("Minutit2","Migrad"))
    results.Print()

    #bkgRange = ROOT.RooRealVar("bkgRange","bkgRange",5.25,5.45) ;
    #bkgRangeArgSet = ROOT.RooArgSet(bkgRange)
    x.setRange("window",5.25,5.45) ;
    fracBkgRange = bkg.createIntegral(obs,obs,"window") ;

    #fracBkgRange = bkg.createIntegral(bkgRangeArgSet,"window") ;
    nbkgWindow = nbkg.getVal() * fracBkgRange.getVal()
    print(nbkg.getVal(), fracBkgRange.getVal())
    print("Number of signals: %f, Number of background: %f, S/sqrt(S+B): %f"%(nsig.getVal(), nbkgWindow, nsig.getVal()/math.sqrt(nsig.getVal() + nbkgWindow)))




    # Plot results of fit on a different frame
    c2 = ROOT.TCanvas('fig_binnedFit', 'fit', 700, 500)
    c2.SetGrid()
    c2.cd()
    ROOT.gPad.SetLeftMargin(0.10)
    ROOT.gPad.SetRightMargin(0.05)

    xframe = wspace.var('x').frame(RooFit.Title("test"))
    hdata.plotOn(xframe)
    model.plotOn(xframe,RooFit.LineColor(2),RooFit.MoveToBack()) # this will show fit overlay on canvas
    #model.plotOn(xframe)
    #model.paramOn(xframe)
    #model.paramOn(xframe) # this will display the fit parameters on canvas
    #filters.Print("v");
    model.plotOn(xframe, RooFit.Components("bkg"),RooFit.LineStyle(ROOT.kDashed),RooFit.LineColor(ROOT.kMagenta),RooFit.MoveToBack()) ;
    model.plotOn(xframe, RooFit.Components("sig"),RooFit.DrawOption("FL"),RooFit.FillColor(9),RooFit.FillStyle(3004),RooFit.LineStyle(6),RooFit.LineColor(9)) ;
    model.plotOn(xframe,RooFit.VisualizeError(results), RooFit.FillColor(ROOT.kOrange), RooFit.MoveToBack()) # this will show fit overlay on canvas

    xframe.GetXaxis().SetTitleFont(42)
    xframe.GetXaxis().SetLabelFont(42)
    xframe.GetYaxis().SetTitleFont(42)
    xframe.GetYaxis().SetLabelFont(42)
    xframe.GetXaxis().SetTitle("m(K^{+}K^{-}e^{+}e^{-}) [GeV]")
    xframe.SetStats(0)
    xframe.SetMinimum(0)
    xframe.Draw()
    c2.SaveAs('test.png')
    print("="*80)

try:
    main()
except KeyboardInterrupt:
    print()
    print("ciao!")
    print()



