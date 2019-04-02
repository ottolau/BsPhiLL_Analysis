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
parser.add_argument("-p", "--oppcharge", dest="oppcharge", default=False, help="Opposite charge")
args = parser.parse_args()


ROOT.gROOT.SetBatch()

def exec_me(command, dryRun=False):
    print(command)
    if not dryRun:
        os.system(command)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def loadAnalyzer(analyzer):
    ROOT.gROOT.ProcessLine(".L analyzer_class/%s.C+"%(analyzer))
    from ROOT import BsPhiLLTupleTree
    #return BsPhiLLTupleTree(tchain)

def makeTraining(inputArgs):
    ich, fChunk , maxevents = inputArgs
    print("Processing chunk number %i"%(ich))
    tchain = ROOT.TChain(args.ttree)
    for filename in fChunk:
        tchain.Add(filename)
    print('Total number of events: ' + str(tchain.GetEntries()))
    #myAnalyzer = loadAnalyzer(tchain, "MVATrainingSet_muon_nontrigger_background")
    myAnalyzer = BsPhiLLTupleTree(tchain)
    oppCharge = False
    if args.oppcharge: oppCharge = True
    myAnalyzer.Loop(outpath+'/'+args.outputfile+'_subset'+str(ich)+'.root', maxevents, 0.0, oppCharge)


if not args.runparallel:
    tchain = ROOT.TChain(args.ttree)
    with open(args.inputfiles) as filenames:
        for filename in filenames:
            tchain.Add(filename.rstrip('\n'))
            break
    #tchain.Add('ggtree_data.root')
    print('Total number of events: ' + str(tchain.GetEntries()))
    #analyzer = "MVATrainingSet_muon_nontrigger_background"
    analyzer = "MVATrainingSet_electron_background"
    #ROOT.gROOT.ProcessLine(".L analyzer_class/%s.C+"%(analyzer))
    ROOT.gSystem.CompileMacro("analyzer_class/%s.C"%(analyzer)) 
    from ROOT import BsPhiLLTupleTree
    myAnalyzer = BsPhiLLTupleTree(tchain)
    myAnalyzer.Loop(args.outputfile, args.maxevents, 0.0, args.oppcharge)

else:

    outputBase = "/eos/uscms/store/user/klau/BsPhiLL_output"
    outputFolder = "BsPhiEE_MVATraining"
    outpath  = "%s/%s"%(outputBase,outputFolder)
    if not os.path.exists(outpath):
        exec_me("mkdir -p %s"%(outpath), False)

    with open(args.inputfiles) as filenames:
        fileList = [f.rstrip('\n') for f in filenames]
    group   = 200
    # stplie files in to n(group) of chunks
    fChunks= list(chunks(fileList,group))
    print("writing %s jobs for %s"%(len(fChunks),outputFolder))
    enumfChunks = enumerate(fChunks)
    maxEventsPerJobs = int(args.maxevents/len(fChunks))
    print("Events per jobs: %i"%(maxEventsPerJobs))
    inputArgs = [it + (maxEventsPerJobs,)  for it in enumfChunks]

    #analyzer = "MVATrainingSet_muon_nontrigger_background"
    analyzer = "MVATrainingSet_electron_background"
    #ROOT.gROOT.ProcessLine(".L analyzer_class/%s.C+"%(analyzer))
    ROOT.gSystem.CompileMacro("analyzer_class/%s.C"%(analyzer)) 
    from ROOT import BsPhiLLTupleTree
    #loadAnalyzer(analyzer)
    pool = mp.Pool(processes = 6)
    pool.map(makeTraining, inputArgs)
    exec_me("hadd %s/%s %s/%s"%(outpath,args.outputfile+'.root',outpath,args.outputfile+'_subset*.root'))
    exec_me("rm %s/%s"%(outpath,args.outputfile+'_subset*.root'))


