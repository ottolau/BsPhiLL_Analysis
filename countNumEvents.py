import ROOT
import os
import multiprocessing as mp

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="DoubleMuonNtu_Run2016B.list", help="List of input ggNtuplizer files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-m", "--maxevents", dest="maxevents", type=int, default=ROOT.TTree.kMaxEntries, help="Maximum number events to loop over")
parser.add_argument("-t", "--ttree", dest="ttree", default="ggNtuplizer/EventTree", help="TTree Name")
parser.add_argument("-r", "--runparallel", dest="runparallel", default=False, help="Enable PROOF")
parser.add_argument("-s", "--sandbox", dest="sandbox", default="", help="Sand Box")
args = parser.parse_args()


ROOT.gROOT.SetBatch()


tchain = ROOT.TChain(args.ttree)
with open(args.inputfiles) as filenames:
    for filename in filenames:
        tchain.Add(filename.rstrip('\n'))
#tchain.Add('ggtree_data.root')
print('Total number of events: ' + str(tchain.GetEntries()))


