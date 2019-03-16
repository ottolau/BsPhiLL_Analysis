#!/usr/bin/env python

import ROOT
from ROOT import TMVA, TFile, TTree, TCut
#from ROOT.TMVA import MethodRXGB
#from ROOT.TMVA import RMethodBase
from subprocess import call
from os.path import isfile

from keras.models import Sequential
from keras.layers.core import Dense, Activation, Dropout
from keras.regularizers import l2
from keras import initializers
from keras.optimizers import SGD, Adam
from keras.layers.normalization import BatchNormalization

#ROOT.gROOT.SetBatch(True)

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
#r = ROOT.TRInterface.Instance()
#TMVA.RMethodBase.Initialize()
#MethodRXGB.Init()

output = TFile.Open('Output_Classification_BsPhiEE.root', 'RECREATE')
#factory = TMVA.Factory('TMVAClassification_BsPhiEE', output,
#        '!V:ROC:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=Classification')

factory = TMVA.Factory('TMVAClassification_BsPhiEE', output,
        '!V:ROC:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

bkg_name = "/eos/uscms/store/user/klau/BsPhiLL_output/BsPhiEE_MVATraining/BsPhiEE_MVATraining_Bkg.root"
sig_name = "/eos/uscms/store/user/klau/BsPhiLL_output/BsPhiEE_MVATraining/BsPhiEE_MVATraining_Sig.root"

branches = ['elePtLead', 'elePtSublead', 'kaonPtLead', 'kaonPtSublead', 'jpsiPt', 'phiPt', 'bsPt', 'eledR', 'kaondR', 'jpsiPhidR', 'svProb', 'svCosine', 'svLxySig', 'eleD0Lead', 'eleD0Sublead', 'eleDzLead', 'eleDzSublead', 'kaonD0Lead', 'kaonD0Sublead', 'kaonDzLead', 'kaonDzSublead', 'kaonNormChi2Lead', 'kaonNormChi2Sublead']

input_dim = len(branches)

use = {}
use['BDT'] = True
use['PyGTB'] = True
use['PyRandomForest'] = True
use['PyAdaBoost'] = True
use['RXGB'] = False
use['PyKeras']  = True
use['DL_DENSE'] = False

data_bkg = TFile.Open(bkg_name)
data_sig = TFile.Open(sig_name)

signal = data_sig.Get('signal')
background = data_bkg.Get('background')

dataloader = TMVA.DataLoader('dataset')

for branch in branches:
    dataloader.AddVariable(branch)

dataloader.AddSignalTree(signal, 1.0)
dataloader.AddBackgroundTree(background, 1.0)
dataloader.PrepareTrainingAndTestTree(TCut(''),
        'nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V')

# Keras DNN

# Define initialization
def normal(shape, name=None):
    return initializers.normal(shape, scale=0.05, name=name)

# Define model
model = Sequential()
model.add(Dense(units=256, input_shape=(input_dim,), kernel_initializer='glorot_normal', activation='relu', kernel_regularizer=l2(1e-5)))
#model.add(Dropout(0.2))
model.add(Dense(units=128, input_shape=(input_dim,), kernel_initializer='glorot_normal', activation='relu', kernel_regularizer=l2(1e-5)))
#model.add(Dropout(0.2))
model.add(Dense(units=64, input_shape=(input_dim,), kernel_initializer='glorot_normal', activation='relu', kernel_regularizer=l2(1e-5)))
#model.add(Dropout(0.2))
model.add(Dense(units=32, input_shape=(input_dim,), kernel_initializer='glorot_normal', activation='relu', kernel_regularizer=l2(1e-5)))
#model.add(Dropout(0.2))
model.add(Dense(2, activation='softmax'))


# Set loss and optimizer
model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=1e-3), metrics=['accuracy',])

# Store model to file
model.save('model.h5')
model.summary()

#TMVA DNN

inputLayoutString = "InputLayout=1|1|%s"%(str(input_dim))
batchLayoutString = "BatchLayout=1|16|%s"%(str(input_dim))

layoutString = "Layout=DENSE|16|RELU,DENSE|16|RELU,DENSE|16|RELU,DENSE|16|RELU,DENSE|1|LINEAR"

training1 = "LearningRate=1e-2,Momentum=0.9,Repetitions=1,ConvergenceSteps=20,BatchSize=16,TestRepetitions=1,MaxEpochs=10,WeightDecay=1e-4,Regularization=L2,Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0."

training2 = "LearningRate=1e-3,Momentum=0.9,Repetitions=1,ConvergenceSteps=20,BatchSize=16,TestRepetitions=1,MaxEpochs=10,WeightDecay=1e-4,Regularization=L2,Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0."

trainingStrategyString = "TrainingStrategy=" + training1 + "|" + training2

dnnOptions = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:WeightInitialization=XAVIERUNIFORM"
dnnOptions = dnnOptions + ":" + inputLayoutString + ":" + batchLayoutString + ":" + layoutString + ":" + trainingStrategyString + ":Architecture=Standard"

#print(dnnOptions)

# Book methods

## TMVA Boosted Decision Trees
if use['BDT']:
    factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDT',
            '!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20')

## Scikit-learn Gradient-Boosted Trees
if use['PyGTB']:
    factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyGTB, "PyGTB","H:V:VarTransform=None:NEstimators=800:LearningRate=0.1:"
                                                      "MaxDepth=3")
## Scikit-learn Random Forest
if use['PyRandomForest']:
    factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyRandomForest, "PyRandomForest","H:V:VarTransform=None:NEstimators=800:"
                               "Criterion=gini:MaxFeatures=auto:MaxDepth=6:MinSamplesLeaf=3:MinWeightFractionLeaf=0:"
                                "Bootstrap=kTRUE" )

## Scikit-learn AdaBoosted Trees      
if use['PyAdaBoost']:
    factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyAdaBoost, "PyAdaBoost","H:V:VarTransform=None:NEstimators=800" )

## XGBoost
if use['RXGB']:
    factory.BookMethod(dataloader, ROOT.TMVA.Types.kRXGB, "RXGB", "!V:NRounds=80:MaxDepth=2:Eta=1" )

## Keras Neural Network
if use['PyKeras']:
    factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras',
            'H:!V:FilenameModel=model.h5:NumEpochs=200:BatchSize=16:VarTransform=None')

## TMVA Neural Network
if use['DL_DENSE']:
    factory.BookMethod(dataloader, TMVA.Types.kDL, 'DL_DENSE', dnnOptions)


# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
