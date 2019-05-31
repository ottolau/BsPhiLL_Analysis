from os import listdir
from os.path import isfile, join

'''
mypath1 = ['/eos/uscms/store/user/klau/ParkingBPH1/crab_ParkingBPH1_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_165252/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH2/crab_ParkingBPH2_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_165826/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH3/crab_ParkingBPH3_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_165914/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH4/crab_ParkingBPH4_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_170014/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH5/crab_ParkingBPH5_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_170102/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH6/crab_ParkingBPH6_Run2018A-14May2018-v1_AOD_18Feb19_BsPhiLL_newTree/190218_170141/0000/']

mypath2 = ['/eos/uscms/store/user/klau/ParkingBPH1/crab_ParkingBPH1_Run2018B-PromptReco-v1_AOD_02Mar19_BsPhiLL_newTree/190303_033546/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH2/crab_ParkingBPH2_Run2018B-PromptReco-v1_AOD_02Mar19_BsPhiLL_newTree/190303_033627/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH3/crab_ParkingBPH3_Run2018B-PromptReco-v1_AOD_02Mar19_BsPhiLL_newTree/190303_033802/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH4/crab_ParkingBPH4_Run2018B-PromptReco-v1_AOD_02Mar19_BsPhiLL_newTree/190303_033845/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH5/crab_ParkingBPH5_Run2018B-PromptReco-v1_AOD_02Mar19_BsPhiLL_newTree/190303_033918/0000/']

mypath3 = ['/eos/uscms/store/user/klau/ParkingBPH1/crab_ParkingBPH1_Run2018D-PromptReco-v2_AOD_03Mar19_BsPhiLL_newTree/190304_041751/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH2/crab_ParkingBPH2_Run2018D-PromptReco-v2_AOD_03Mar19_BsPhiLL_newTree/190304_041941/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH3/crab_ParkingBPH3_Run2018D-PromptReco-v2_AOD_03Mar19_BsPhiLL_newTree/190304_042301/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH4/crab_ParkingBPH4_Run2018D-PromptReco-v2_AOD_03Mar19_BsPhiLL_newTree/190304_042411/0000/',
           '/eos/uscms/store/user/klau/ParkingBPH5/crab_ParkingBPH5_Run2018D-PromptReco-v2_AOD_03Mar19_BsPhiLL_newTree/190304_042529/0000/']

mypath = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_02Apr19_BsPhiLL_lowPtElectrons/190403_014932/0000/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_02Apr19_BsPhiLL_lowPtElectrons/190403_014932/0001/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_02Apr19_BsPhiLL_lowPtElectrons/190403_014932/0002/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_02Apr19_BsPhiLL_lowPtElectrons/190403_014932/0003/']

mypath1 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190405_003242/0000/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190405_003242/0001/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190405_003242/0002/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190405_003242/0003/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190405_003242/0004/']


mypath2 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190404_191546/0000/',
          '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_04Apr19_BsPhiLL_lowPtElectrons/190404_191546/0001/']


mypath1 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181324/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181324/0001/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181324/0002/']

mypath2 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018D-20Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181412/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH2/crab_ParkingBPH2_Run2018D-20Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181631/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181709/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH4/crab_ParkingBPH4_Run2018D-20Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181822/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH5/crab_ParkingBPH5_Run2018D-20Mar2019-v1_AOD_10Apr19_BsPhiLL_lowPtElectrons_ptMode/190410_181959/0000/']


mypath1 = ['/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_184708/0000/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_184708/0001/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_184708/0002/']

mypath2 = ['/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH1/crab_ParkingBPH1_Run2018D-20Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_185905/0000/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH2/crab_ParkingBPH2_Run2018D-20Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_190008/0000/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_190116/0000/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH4/crab_ParkingBPH4_Run2018D-20Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_190144/0000/',
           '/eos/uscms/store/user/klau/DielectronSpectrum/ParkingBPH5/crab_ParkingBPH5_Run2018D-20Mar2019-v1_AOD_10Apr19_DielectronSpectrum_ptMode/190410_190240/0000/']


mypath1 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210044/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210044/0001/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210044/0002/']

mypath2 = ['/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018D-20Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210508/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH2/crab_ParkingBPH2_Run2018D-20Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210610/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210711/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH4/crab_ParkingBPH4_Run2018D-20Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210750/0000/',
           '/eos/uscms/store/user/klau/LowPtElectrons/ParkingBPH5/crab_ParkingBPH5_Run2018D-20Mar2019-v1_AOD_14Apr19_BsPhiLL_lowPtElectrons_UnBWPVeryLoose/190414_210828/0000/']

mypath1 = ['/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233013/0000/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233013/0001/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233013/0002/']

mypath2 = ['/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH1/crab_ParkingBPH1_Run2018D-20Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233318/0000/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH2/crab_ParkingBPH2_Run2018D-20Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233646/0000/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233717/0000/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH4/crab_ParkingBPH4_Run2018D-20Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233801/0000/',
           '/eos/uscms/store/user/klau/PFAfterSkim/ParkingBPH5/crab_ParkingBPH5_Run2018D-20Mar2019-v1_AOD_15Apr19_PFAfterSkim/190415_233852/0000/']
'''
mypath1 = ['/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_190757/0000/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_190757/0001/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018A-22Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_190757/0002/']

mypath2 = ['/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH1/crab_ParkingBPH1_Run2018D-20Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_192147/0000/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH2/crab_ParkingBPH2_Run2018D-20Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_192220/0000/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH3/crab_ParkingBPH3_Run2018D-20Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_192334/0000/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH4/crab_ParkingBPH4_Run2018D-20Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_192447/0000/',
           '/eos/uscms/store/user/klau/BKLL_LowPtElectrons/ParkingBPH5/crab_ParkingBPH5_Run2018D-20Mar2019-v1_AOD_16Apr19_BKLL_lowPtElectrons/190416_192515/0000/']


mypath = mypath1 + mypath2 #+ mypath3


filename = []
for path in mypath:
    filename = filename + [path + f for f in listdir(path) if isfile(join(path, f))]
redirector = 'root://cmsxrootd.fnal.gov/'
filename = ['%s%s'%(redirector, f[10:]) for f in filename if ".root" in f]

#outputFile = open('ParkingBPH123456Ntu_Run2018A_4tracks_MINIAOD_newTree.list', 'w+')
#outputFile = open('ParkingBPH12345Ntu_Run2018B_4tracks_MINIAOD_newTree.list', 'w+')
#outputFile = open('ParkingBPH2345Ntu_Run2018D_4tracks_MINIAOD_newTree.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018ABD_4tracks_MINIAOD_newTree.list', 'w+')
#outputFile = open('ParkingBPH123456Ntu_Run2018A_AOD_BsPhiLL_newTree.list', 'w+')
#outputFile = open('ParkingBPH12345Ntu_Run2018B_AOD_BsPhiLL_newTree.list', 'w+')
#outputFile = open('ParkingBPH12345Ntu_Run2018D_AOD_BsPhiLL_newTree.list', 'w+')
#outputFile = open('ParkingBPH1Ntu_Run2018A_AOD_BsPhiLL_lowPtElectrons.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D3_AOD_BsPhiLL_lowPtElectrons.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D3_AOD_BsPhiLL_dilectron.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D12345_AOD_BsPhiLL_lowPtElectron_ptMode.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D12345_AOD_BsPhiLL_dilectron_ptMode.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D12345_AOD_BsPhiLL_lowPtElectron_UnBWPVeryLoose.list', 'w+')
#outputFile = open('ParkingBPHNtu_Run2018A1D12345_AOD_BsPhiLL_PFAfterSkim.list', 'w+')
outputFile = open('ParkingBPHNtu_Run2018A1D12345_AOD_BKLL_lowPtElectron.list', 'w+')


for f in filename:
    outputFile.write('%s\n'%(f))
outputFile.close()

