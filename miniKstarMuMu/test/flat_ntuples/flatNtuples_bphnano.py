from argparse import ArgumentParser
from samples import *
import pdb 

parser = ArgumentParser()
# parser.add_argument("sample",                                       help = "sample", nargs = '+', choices = samples.keys(), default = 'BJpsiK_ee_mc_2019Oct25' )
parser.add_argument("-n"  , "--number"     , dest = "number"     ,  help = "number of file"              , default = '1')
parser.add_argument("-f"  , "--folder"     , dest = "folder"     ,  help = "name of output folder"       , default = 'default_folder')
# parser.add_argument("-d"  , "--dir"        , dest = "indir"      ,  help = "name of input folder (0000)" , default = '0000')
# parser.add_argument("-e"  , "--era"        , dest = "theera"     ,  help = "2018B,2018B..."              , default = '2016B')
# parser.add_argument("-c"  , "--chan"       , dest = "thechan"    ,  help = "LMNR, PSI"                   , default = 'LMNR' )

args = parser.parse_args()
if not args.number:   
  parser.error('Number filename not given')
  
import ROOT
ROOT.gSystem.Load('libminiB0KstarMuMuminiKstarMuMu')
# from ROOT               import miniHLTObj
from DataFormats.FWLite import Events
from collections        import OrderedDict
from array              import array
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

import sys, math, itertools
sys.path.append('/gwpool/users/fiorendi/p5prime/run3/CMSSW_14_2_2/src/miniB0KstarMuMu/miniKstarMuMu')
from utils.angular_vars import *
from utils.utils import *

import numpy as np
from ROOT import TLorentzVector

###############################################
## Parameters for the flattening:

# type       = args.thechan 
skim       = True
skimSoftMu = True
# isNot2016  = False if '2016' in args.sample[0] else True
paths = [ 
          'HLT_DoubleMu4_LowMassNonResonantTrk_Displaced',
        ]  

          
###############################################
tree_lmnr = ROOT.TChain('Events')

tree_lmnr.Add('/gwpool/users/fiorendi/p5prime/run3/CMSSW_14_2_2/src/miniB0KstarMuMu/miniKstarMuMu/test/flat_ntuples/BPH_newversion.root')
file_out = ROOT.TFile( 'ntuple_flat.root', 'recreate')
ntuple   = ROOT.TTree( 'ntuple', 'ntuple' )


print ('@@@@@@@@@@@     flattening ntuple     @@@@@@@@@@@')
print ('input file: ', tree_lmnr.GetFile())

try:    tree_lmnr.GetEntries()
except:
  print ('ntuple tree not found')
  exit()


newbr = []
isobr = []
evbr  = []
# 
# '''
# for i in branches:
#     print "{NAME}\t = array('d', [-99.]);  newbr.append({NAME})".format(NAME=i).expandtabs(40)
# '''
eventN                           = array('L', [ 99 ]);  
runN                             = array('d', [-99.]);  evbr .append(runN)
recoVtxN                         = array('d', [-99.]);  evbr .append(recoVtxN)
bsX                              = array('d', [-99.]);  evbr .append(bsX)
bsY                              = array('d', [-99.]);  evbr .append(bsY)

# trig                             = array('d', [-99.]);  newbr.append(trig)
# l1_00_1p5                        = array('d', [-99.]);  newbr.append(l1_00_1p5 )
# l1_00_1p4                        = array('d', [-99.]);  newbr.append(l1_00_1p4 )
# l1_4_4                           = array('d', [-99.]);  newbr.append(l1_4_4    )
# l1_4p5_4p5                       = array('d', [-99.]);  newbr.append(l1_4p5_4p5)
# # l1_11_4                          = array('d', [-99.]);  newbr.append(l1_11_4 )
# # l1_12_5                          = array('d', [-99.]);  newbr.append(l1_12_5 )
# # l1_10_0                          = array('d', [-99.]);  newbr.append(l1_10_0 )
# # l1_00                            = array('d', [-99.]);  newbr.append(l1_00   )
# # l1_00_OS                         = array('d', [-99.]);  newbr.append(l1_00_OS)
bMass                            = array('d', [-99.]);  newbr.append(bMass)
bMassE                           = array('d', [-99.]);  newbr.append(bMassE)
bBarMass                         = array('d', [-99.]);  newbr.append(bBarMass)
bBarMassE                        = array('d', [-99.]);  newbr.append(bBarMassE)
bVtxCL                           = array('d', [-99.]);  newbr.append(bVtxCL)
bVtxX                            = array('d', [-99.]);  newbr.append(bVtxX)
bVtxY                            = array('d', [-99.]);  newbr.append(bVtxY)
bVtxZ                            = array('d', [-99.]);  newbr.append(bVtxZ)
bCosAlphaBS                      = array('d', [-99.]);  newbr.append(bCosAlphaBS)
# bCosAlphaBSE                     = array('d', [-99.]);  newbr.append(bCosAlphaBSE)
bLBS                             = array('d', [-99.]);  newbr.append(bLBS)
bLBSE                            = array('d', [-99.]);  newbr.append(bLBSE)
# bDCABS                           = array('d', [-99.]);  newbr.append(bDCABS)
# bDCABSE                          = array('d', [-99.]);  newbr.append(bDCABSE)
# 
kstMass                          = array('d', [-99.]);  newbr.append(kstMass)
kstMassE                         = array('d', [-99.]);  newbr.append(kstMassE)
kstBarMass                       = array('d', [-99.]);  newbr.append(kstBarMass)
kstBarMassE                      = array('d', [-99.]);  newbr.append(kstBarMassE)
kstVtxCL                         = array('d', [-99.]);  newbr.append(kstVtxCL)
kstVtxX                          = array('d', [-99.]);  newbr.append(kstVtxX)
kstVtxY                          = array('d', [-99.]);  newbr.append(kstVtxY)
kstVtxZ                          = array('d', [-99.]);  newbr.append(kstVtxZ)
kkMass                           = array('d', [-99.]);  newbr.append(kkMass)

mumuMass                         = array('d', [-99.]);  newbr.append(mumuMass)
# mumuMassE                        = array('d', [-99.]);  newbr.append(mumuMassE)
# mumuVtxCL                        = array('d', [-99.]);  newbr.append(mumuVtxCL)
# mumuVtxX                         = array('d', [-99.]);  newbr.append(mumuVtxX)
# mumuVtxY                         = array('d', [-99.]);  newbr.append(mumuVtxY)
# mumuVtxZ                         = array('d', [-99.]);  newbr.append(mumuVtxZ)
# mumuCosAlphaBS                   = array('d', [-99.]);  newbr.append(mumuCosAlphaBS)
# mumuCosAlphaBSE                  = array('d', [-99.]);  newbr.append(mumuCosAlphaBSE)
# mumuLBS                          = array('d', [-99.]);  newbr.append(mumuLBS)
# mumuLBSE                         = array('d', [-99.]);  newbr.append(mumuLBSE)
# mumuDCA                          = array('d', [-99.]);  newbr.append(mumuDCA)
# 

mumSoft                = array('d', [-99.]);  newbr.append(mumSoft)
mumSoftMVA             = array('d', [-99.]);  newbr.append(mumSoftMVA)
mumSoftMVARun3Score    = array('d', [-99.]);  newbr.append(mumSoftMVARun3Score)

mupSoft                = array('d', [-99.]);  newbr.append(mupSoft)
mupSoftMVA             = array('d', [-99.]);  newbr.append(mupSoftMVA)
mupSoftMVARun3Score    = array('d', [-99.]);  newbr.append(mupSoftMVARun3Score)


# mumHighPurity                    = array('d', [-99.]);  newbr.append(mumHighPurity)
# mumCL                            = array('d', [-99.]);  newbr.append(mumCL)
# mumNormChi2                      = array('d', [-99.]);  newbr.append(mumNormChi2)
# mumDCAVtx                        = array('d', [-99.]);  newbr.append(mumDCAVtx)
# mumDCAVtxE                       = array('d', [-99.]);  newbr.append(mumDCAVtxE)
# mumdxyBS                         = array('d', [-99.]);  newbr.append(mumdxyBS)
# mumdzBS                          = array('d', [-99.]);  newbr.append(mumdzBS)
# mumMinIP2D                       = array('d', [-99.]);  newbr.append(mumMinIP2D)
# mumMinIP2DE                      = array('d', [-99.]);  newbr.append(mumMinIP2DE)
# mumNPixLayers                    = array('d', [-99.]);  newbr.append(mumNPixLayers)
# mumNTrkLayers                    = array('d', [-99.]);  newbr.append(mumNTrkLayers)
# mumNMuonHits                     = array('d', [-99.]);  newbr.append(mumNMuonHits)
# mumNMatchStation                 = array('d', [-99.]);  newbr.append(mumNMatchStation)
# mupHighPurity                    = array('d', [-99.]);  newbr.append(mupHighPurity)
# mupCL                            = array('d', [-99.]);  newbr.append(mupCL)
# mupNormChi2                      = array('d', [-99.]);  newbr.append(mupNormChi2)
# mupDCAVtx                        = array('d', [-99.]);  newbr.append(mupDCAVtx)
# mupDCAVtxE                       = array('d', [-99.]);  newbr.append(mupDCAVtxE)
# mupdxyBS                         = array('d', [-99.]);  newbr.append(mupdxyBS)
# mupdzBS                          = array('d', [-99.]);  newbr.append(mupdzBS)
# mupMinIP2D                       = array('d', [-99.]);  newbr.append(mupMinIP2D)
# mupMinIP2DE                      = array('d', [-99.]);  newbr.append(mupMinIP2DE)
# mupNPixLayers                    = array('d', [-99.]);  newbr.append(mupNPixLayers)
# mupNTrkLayers                    = array('d', [-99.]);  newbr.append(mupNTrkLayers)
# mupNMuonHits                     = array('d', [-99.]);  newbr.append(mupNMuonHits)
# mupNMatchStation                 = array('d', [-99.]);  newbr.append(mupNMatchStation)
# kstTrkmHighPurity                = array('d', [-99.]);  newbr.append(kstTrkmHighPurity)
# kstTrkmCL                        = array('d', [-99.]);  newbr.append(kstTrkmCL)
# kstTrkmNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkmNormChi2)
# kstTrkmPz                        = array('d', [-99.]);  newbr.append(kstTrkmPz)
kstTrkmDCABSSign                     = array('d', [-99.]);  newbr.append(kstTrkmDCABSSign)
# kstTrkmDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkmDCABSE)
# kstTrkmFracHits                  = array('d', [-99.]);  newbr.append(kstTrkmFracHits)
# kstTrkmMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkmMinIP2D)
# kstTrkmMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkmMinIP2DE)
kstTrkmNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkmNPixHits)
# kstTrkmNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkmNPixLayers)
# kstTrkmNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkmNTrkHits)
# kstTrkmNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkmNTrkLayers)
# kstTrkmMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkmMuMatch)
# kstTrkpHighPurity                = array('d', [-99.]);  newbr.append(kstTrkpHighPurity)
# kstTrkpCL                        = array('d', [-99.]);  newbr.append(kstTrkpCL)
# kstTrkpNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkpNormChi2)
# kstTrkpPz                        = array('d', [-99.]);  newbr.append(kstTrkpPz)
kstTrkpDCABSSign                     = array('d', [-99.]);  newbr.append(kstTrkpDCABSSign)
# kstTrkpDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkpDCABSE)
# kstTrkpFracHits                  = array('d', [-99.]);  newbr.append(kstTrkpFracHits)
# kstTrkpMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkpMinIP2D)
# kstTrkpMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkpMinIP2DE)
kstTrkpNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkpNPixHits)
# kstTrkpNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkpNPixLayers)
# kstTrkpNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkpNTrkHits)
# kstTrkpNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkpNTrkLayers)
# kstTrkpMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkpMuMatch)
tagB0                            = array('d', [-99.]);  newbr.append(tagB0)
# 
# bMinusVtxCL                      = array('d', [-99.]);  newbr.append(bMinusVtxCL      )
# bMinusCosAlphaBS                 = array('d', [-99.]);  newbr.append(bMinusCosAlphaBS )
# bPlusVtxCL                       = array('d', [-99.]);  newbr.append(bPlusVtxCL       )
# bPlusCosAlphaBS                  = array('d', [-99.]);  newbr.append(bPlusCosAlphaBS  )
# 
mumGlobalMuon                    = array('d',[-99.]); newbr.append(mumGlobalMuon             )   
mumTrackerMuon                   = array('d',[-99.]); newbr.append(mumTrackerMuon            )   
# mumStandAloneMuon                = array('d',[-99.]); newbr.append(mumStandAloneMuon         )   
# mumTMOneStationTight             = array('d',[-99.]); newbr.append(mumTMOneStationTight      )   
# mumTMOneStationLoose             = array('d',[-99.]); newbr.append(mumTMOneStationLoose      )   
mupGlobalMuon                    = array('d',[-99.]); newbr.append(mupGlobalMuon             )   
mupTrackerMuon                   = array('d',[-99.]); newbr.append(mupTrackerMuon            )   
# mupStandAloneMuon                = array('d',[-99.]); newbr.append(mupStandAloneMuon         )   
# mupTMOneStationTight             = array('d',[-99.]); newbr.append(mupTMOneStationTight      )   
# mupTMOneStationLoose             = array('d',[-99.]); newbr.append(mupTMOneStationLoose      )   
# kstTrkmGlobalMuon                = array('d',[-99.]); newbr.append(kstTrkmGlobalMuon             )   
# kstTrkmTrackerMuon               = array('d',[-99.]); newbr.append(kstTrkmTrackerMuon            )   
# kstTrkmStandAloneMuon            = array('d',[-99.]); newbr.append(kstTrkmStandAloneMuon         )   
# kstTrkmGlobalMuonPromptTight     = array('d',[-99.]); newbr.append(kstTrkmGlobalMuonPromptTight  )   
# kstTrkmTMOneStationTight         = array('d',[-99.]); newbr.append(kstTrkmTMOneStationTight      )   
# kstTrkmTMOneStationLoose         = array('d',[-99.]); newbr.append(kstTrkmTMOneStationLoose      )   
# kstTrkmTrackerMuonArbitrated     = array('d',[-99.]); newbr.append(kstTrkmTrackerMuonArbitrated  )   
# kstTrkpGlobalMuon                = array('d',[-99.]); newbr.append(kstTrkpGlobalMuon             )   
# kstTrkpTrackerMuon               = array('d',[-99.]); newbr.append(kstTrkpTrackerMuon            )   
# kstTrkpStandAloneMuon            = array('d',[-99.]); newbr.append(kstTrkpStandAloneMuon         )   
# kstTrkpGlobalMuonPromptTight     = array('d',[-99.]); newbr.append(kstTrkpGlobalMuonPromptTight  )   
# kstTrkpTMOneStationTight         = array('d',[-99.]); newbr.append(kstTrkpTMOneStationTight      )   
# kstTrkpTMOneStationLoose         = array('d',[-99.]); newbr.append(kstTrkpTMOneStationLoose      )   
# kstTrkpTrackerMuonArbitrated     = array('d',[-99.]); newbr.append(kstTrkpTrackerMuonArbitrated  )   
# 
# 
bPt                              = array('d',[-99.]); newbr.append( bPt          )   
kstPt                            = array('d',[-99.]); newbr.append( kstPt        )   
mumuPt                           = array('d',[-99.]); newbr.append( mumuPt       )   
mumPt                            = array('d',[-99.]); newbr.append( mumPt        )   
mupPt                            = array('d',[-99.]); newbr.append( mupPt        )   
mupCharge                        = array('d',[-99.]); newbr.append( mupCharge        )   ## test for BPHnano
kstTrkmPt                        = array('d',[-99.]); newbr.append( kstTrkmPt    )   
kstTrkpPt                        = array('d',[-99.]); newbr.append( kstTrkpPt    )   
bPhi                             = array('d',[-99.]); newbr.append( bPhi         )   
kstPhi                           = array('d',[-99.]); newbr.append( kstPhi       )   
mumuPhi                          = array('d',[-99.]); newbr.append( mumuPhi      )   
mumPhi                           = array('d',[-99.]); newbr.append( mumPhi       )   
mupPhi                           = array('d',[-99.]); newbr.append( mupPhi       )   
kstTrkmPhi                       = array('d',[-99.]); newbr.append( kstTrkmPhi   )   
kstTrkpPhi                       = array('d',[-99.]); newbr.append( kstTrkpPhi   )   
bEta                             = array('d',[-99.]); newbr.append( bEta         )   
kstEta                           = array('d',[-99.]); newbr.append( kstEta       )   
mumuEta                          = array('d',[-99.]); newbr.append( mumuEta      )   
mumEta                           = array('d',[-99.]); newbr.append( mumEta       )   
mupEta                           = array('d',[-99.]); newbr.append( mupEta       )   
kstTrkmEta                       = array('d',[-99.]); newbr.append( kstTrkmEta   )   
kstTrkpEta                       = array('d',[-99.]); newbr.append( kstTrkpEta   )   
# 
# mmk1                              = array('d',[-99.]); newbr.append( mmk1          )   
# mmk2                              = array('d',[-99.]); newbr.append( mmk2          )   
dR_mum_trkm                       = array('d',[-99.]); newbr.append( dR_mum_trkm   )          
dR_mup_trkp                       = array('d',[-99.]); newbr.append( dR_mup_trkp   )          

cos_theta_l                      = array('d', [-99.]); newbr.append(cos_theta_l)
cos_theta_k                      = array('d', [-99.]); newbr.append(cos_theta_k)
phi_kst_mumu                     = array('d', [-99.]); newbr.append(phi_kst_mumu)

# mumIsoPt_dr03                    = array('d',[0.]); isobr.append( mumIsoPt_dr03   )   
# mupIsoPt_dr03                    = array('d',[0.]); isobr.append( mupIsoPt_dr03   )   
# kstTrkmIsoPt_dr03                = array('d',[0.]); isobr.append( kstTrkmIsoPt_dr03   )   
# kstTrkpIsoPt_dr03                = array('d',[0.]); isobr.append( kstTrkpIsoPt_dr03   )   
# 
# mumIsoPt_dr04                    = array('d',[0.]); isobr.append( mumIsoPt_dr04   )   
# mupIsoPt_dr04                    = array('d',[0.]); isobr.append( mupIsoPt_dr04   )   
# kstTrkmIsoPt_dr04                = array('d',[0.]); isobr.append( kstTrkmIsoPt_dr04   )   
# kstTrkpIsoPt_dr04                = array('d',[0.]); isobr.append( kstTrkpIsoPt_dr04   )   
# 
# '''
# for i in branches:
#     print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/D')".format(NAME=i).expandtabs(40)
# '''
ntuple.Branch('runN',                   runN,                                   'runN/D')
ntuple.Branch('eventN',                 eventN,                                 'eventN/L')
ntuple.Branch('recoVtxN',               recoVtxN,                               'recoVtxN/D')
ntuple.Branch('bsX',                    bsX,                                    'bsX/D')
ntuple.Branch('bsY',                    bsY,                                    'bsY/D')
# 
# ntuple.Branch('trig',                   trig,                                   'trig/D')
# ntuple.Branch('l1_00_1p5',              l1_00_1p5 ,                             'l1_00_1p5/D')
# ntuple.Branch('l1_00_1p4',              l1_00_1p4 ,                             'l1_00_1p4/D')
# ntuple.Branch('l1_4_4',                 l1_4_4 ,                                'l1_4_4/D')
# ntuple.Branch('l1_4p5_4p5',             l1_4p5_4p5 ,                            'l1_4p5_4p5/D')
# 
ntuple.Branch('bMass',                  bMass,                                  'bMass/D')
ntuple.Branch('bMassE',                 bMassE,                                 'bMassE/D')
ntuple.Branch('bBarMass',               bBarMass,                               'bBarMass/D')
ntuple.Branch('bBarMassE',              bBarMassE,                              'bBarMassE/D')
ntuple.Branch('bVtxCL',                 bVtxCL,                                 'bVtxCL/D')
ntuple.Branch('bVtxX',                  bVtxX,                                  'bVtxX/D')
ntuple.Branch('bVtxY',                  bVtxY,                                  'bVtxY/D')
ntuple.Branch('bVtxZ',                  bVtxZ,                                  'bVtxZ/D')
# 
ntuple.Branch('bCosAlphaBS',            bCosAlphaBS,                            'bCosAlphaBS/D')
# ntuple.Branch('bCosAlphaBSE',           bCosAlphaBSE,                           'bCosAlphaBSE/D')
ntuple.Branch('bLBS',                   bLBS,                                   'bLBS/D')
ntuple.Branch('bLBSE',                  bLBSE,                                  'bLBSE/D')
# ntuple.Branch('bDCABS',                 bDCABS,                                 'bDCABS/D')
# ntuple.Branch('bDCABSE',                bDCABSE,                                'bDCABSE/D')
# 
ntuple.Branch('kstMass',                kstMass,                                'kstMass/D')
ntuple.Branch('kstMassE',               kstMassE,                               'kstMassE/D')
ntuple.Branch('kstBarMass',             kstBarMass,                             'kstBarMass/D')
ntuple.Branch('kstBarMassE',            kstBarMassE,                            'kstBarMassE/D')
ntuple.Branch('kstVtxCL',               kstVtxCL,                               'kstVtxCL/D')
ntuple.Branch('kstVtxX',                kstVtxX,                                'kstVtxX/D')
ntuple.Branch('kstVtxY',                kstVtxY,                                'kstVtxY/D')
ntuple.Branch('kstVtxZ',                kstVtxZ,                                'kstVtxZ/D')
ntuple.Branch('kkMass',                 kkMass,                                 'kkMass/D')
ntuple.Branch('mumuMass',               mumuMass,                               'mumuMass/D')
# ntuple.Branch('mumuMassE',              mumuMassE,                              'mumuMassE/D')
# ntuple.Branch('mumuVtxCL',              mumuVtxCL,                              'mumuVtxCL/D')
# ntuple.Branch('mumuVtxX',               mumuVtxX,                               'mumuVtxX/D')
# ntuple.Branch('mumuVtxY',               mumuVtxY,                               'mumuVtxY/D')
# ntuple.Branch('mumuVtxZ',               mumuVtxZ,                               'mumuVtxZ/D')
# ntuple.Branch('mumuCosAlphaBS',         mumuCosAlphaBS,                         'mumuCosAlphaBS/D')
# ntuple.Branch('mumuCosAlphaBSE',        mumuCosAlphaBSE,                        'mumuCosAlphaBSE/D')
# ntuple.Branch('mumuLBS',                mumuLBS,                                'mumuLBS/D')
# ntuple.Branch('mumuLBSE',               mumuLBSE,                               'mumuLBSE/D')
# ntuple.Branch('mumuDCA',                mumuDCA,                                'mumuDCA/D')
# 

ntuple.Branch('mumSoft',             mumSoft             ,                          'mumSoft/D')
ntuple.Branch('mumSoftMVA',          mumSoftMVA          ,                          'mumSoftMVA/D')
ntuple.Branch('mumSoftMVARun3Score', mumSoftMVARun3Score ,                          'mumSoftMVARun3Score/D')
ntuple.Branch('mupSoft',             mupSoft             ,                          'mupSoft/D')
ntuple.Branch('mupSoftMVA',          mupSoftMVA          ,                          'mupSoftMVA/D')
ntuple.Branch('mupSoftMVARun3Score', mupSoftMVARun3Score ,                          'mupSoftMVARun3Score/D')

# ntuple.Branch('mumHighPurity',          mumHighPurity,                          'mumHighPurity/D')
# ntuple.Branch('mumCL',                  mumCL,                                  'mumCL/D')
# ntuple.Branch('mumNormChi2',            mumNormChi2,                            'mumNormChi2/D')
# ntuple.Branch('mumDCAVtx',              mumDCAVtx,                              'mumDCAVtx/D')
# ntuple.Branch('mumDCAVtxE',             mumDCAVtxE,                             'mumDCAVtxE/D')
# ntuple.Branch('mumdxyBS',               mumdxyBS,                               'mumdxyBS/D')
# ntuple.Branch('mumdzBS',                mumdzBS,                                'mumdzBS/D')
# ntuple.Branch('mumMinIP2D',             mumMinIP2D,                             'mumMinIP2D/D')
# ntuple.Branch('mumMinIP2DE',            mumMinIP2DE,                            'mumMinIP2DE/D')
# ntuple.Branch('mumNPixLayers',          mumNPixLayers,                          'mumNPixLayers/D')
# ntuple.Branch('mumNTrkLayers',          mumNTrkLayers,                          'mumNTrkLayers/D')
# ntuple.Branch('mumNMuonHits',           mumNMuonHits,                           'mumNMuonHits/D')
# ntuple.Branch('mumNMatchStation',       mumNMatchStation,                       'mumNMatchStation/D')
# ntuple.Branch('mupHighPurity',          mupHighPurity,                          'mupHighPurity/D')
# ntuple.Branch('mupCL',                  mupCL,                                  'mupCL/D')
# ntuple.Branch('mupNormChi2',            mupNormChi2,                            'mupNormChi2/D')
# ntuple.Branch('mupDCAVtx',              mupDCAVtx,                              'mupDCAVtx/D')
# ntuple.Branch('mupDCAVtxE',             mupDCAVtxE,                             'mupDCAVtxE/D')
# ntuple.Branch('mupdxyBS',               mupdxyBS,                               'mupdxyBS/D')
# ntuple.Branch('mupdzBS',                mupdzBS,                                'mupdzBS/D')
# ntuple.Branch('mupMinIP2D',             mupMinIP2D,                             'mupMinIP2D/D')
# ntuple.Branch('mupMinIP2DE',            mupMinIP2DE,                            'mupMinIP2DE/D')
# ntuple.Branch('mupNPixLayers',          mupNPixLayers,                          'mupNPixLayers/D')
# ntuple.Branch('mupNTrkLayers',          mupNTrkLayers,                          'mupNTrkLayers/D')
# ntuple.Branch('mupNMuonHits',           mupNMuonHits,                           'mupNMuonHits/D')
# ntuple.Branch('mupNMatchStation',       mupNMatchStation,                       'mupNMatchStation/D')
# 
# ntuple.Branch('kstTrkmHighPurity',      kstTrkmHighPurity,                      'kstTrkmHighPurity/D')
# ntuple.Branch('kstTrkmCL',              kstTrkmCL,                              'kstTrkmCL/D')
# ntuple.Branch('kstTrkmNormChi2',        kstTrkmNormChi2,                        'kstTrkmNormChi2/D')
# ntuple.Branch('kstTrkmPz',              kstTrkmPz,                              'kstTrkmPz/D')
ntuple.Branch('kstTrkmDCABSSign',           kstTrkmDCABSSign,                           'kstTrkmDCABSSign/D')
# ntuple.Branch('kstTrkmDCABSE',          kstTrkmDCABSE,                          'kstTrkmDCABSE/D')
# ntuple.Branch('kstTrkmFracHits',        kstTrkmFracHits,                        'kstTrkmFracHits/D')
# ntuple.Branch('kstTrkmMinIP2D',         kstTrkmMinIP2D,                         'kstTrkmMinIP2D/D')
# ntuple.Branch('kstTrkmMinIP2DE',        kstTrkmMinIP2DE,                        'kstTrkmMinIP2DE/D')
ntuple.Branch('kstTrkmNPixHits',        kstTrkmNPixHits,                        'kstTrkmNPixHits/D')
# ntuple.Branch('kstTrkmNPixLayers',      kstTrkmNPixLayers,                      'kstTrkmNPixLayers/D')
# ntuple.Branch('kstTrkmNTrkHits',        kstTrkmNTrkHits,                        'kstTrkmNTrkHits/D')
# ntuple.Branch('kstTrkmNTrkLayers',      kstTrkmNTrkLayers,                      'kstTrkmNTrkLayers/D')
# ntuple.Branch('kstTrkmMuMatch',         kstTrkmMuMatch,                         'kstTrkmMuMatch/D')
# ntuple.Branch('kstTrkpHighPurity',      kstTrkpHighPurity,                      'kstTrkpHighPurity/D')
# ntuple.Branch('kstTrkpCL',              kstTrkpCL,                              'kstTrkpCL/D')
# ntuple.Branch('kstTrkpNormChi2',        kstTrkpNormChi2,                        'kstTrkpNormChi2/D')
# ntuple.Branch('kstTrkpPz',              kstTrkpPz,                              'kstTrkpPz/D')
ntuple.Branch('kstTrkpDCABSSign',           kstTrkpDCABSSign,                           'kstTrkpDCABSSign/D')
# ntuple.Branch('kstTrkpDCABSE',          kstTrkpDCABSE,                          'kstTrkpDCABSE/D')
# ntuple.Branch('kstTrkpFracHits',        kstTrkpFracHits,                        'kstTrkpFracHits/D')
# ntuple.Branch('kstTrkpMinIP2D',         kstTrkpMinIP2D,                         'kstTrkpMinIP2D/D')
# ntuple.Branch('kstTrkpMinIP2DE',        kstTrkpMinIP2DE,                        'kstTrkpMinIP2DE/D')
ntuple.Branch('kstTrkpNPixHits',        kstTrkpNPixHits,                        'kstTrkpNPixHits/D')
# ntuple.Branch('kstTrkpNPixLayers',      kstTrkpNPixLayers,                      'kstTrkpNPixLayers/D')
# ntuple.Branch('kstTrkpNTrkHits',        kstTrkpNTrkHits,                        'kstTrkpNTrkHits/D')
# ntuple.Branch('kstTrkpNTrkLayers',      kstTrkpNTrkLayers,                      'kstTrkpNTrkLayers/D')
# ntuple.Branch('kstTrkpMuMatch',         kstTrkpMuMatch,                         'kstTrkpMuMatch/D')
ntuple.Branch('tagB0',                  tagB0,                                  'tagB0/D')
ntuple.Branch('cos_theta_l' ,           cos_theta_l ,                           'cos_theta_l/D')
ntuple.Branch('cos_theta_k' ,           cos_theta_k ,                           'cos_theta_k/D')
ntuple.Branch('phi_kst_mumu',           phi_kst_mumu,                           'phi_kst_mumu/D')
# 
# ntuple.Branch('bMinusVtxCL',            bMinusVtxCL     ,                       'bMinusVtxCL/D')
# ntuple.Branch('bMinusCosAlphaBS',       bMinusCosAlphaBS,                       'bMinusCosAlphaBS/D')
# ntuple.Branch('bPlusVtxCL',             bPlusVtxCL      ,                       'bPlusVtxCL/D')
# ntuple.Branch('bPlusCosAlphaBS',        bPlusCosAlphaBS ,                       'bPlusCosAlphaBS/D')
# 
ntuple.Branch('mumGlobalMuon',             mumGlobalMuon              ,               'mumGlobalMuon/D')
ntuple.Branch('mumTrackerMuon',            mumTrackerMuon             ,               'mumTrackerMuon/D')
# ntuple.Branch('mumStandAloneMuon',         mumStandAloneMuon          ,               'mumStandAloneMuon/D')
# ntuple.Branch('mumTMOneStationTight',      mumTMOneStationTight       ,               'mumTMOneStationTight/D')
# ntuple.Branch('mumTMOneStationLoose',      mumTMOneStationLoose       ,               'mumTMOneStationLoose/D')
ntuple.Branch('mupGlobalMuon',             mupGlobalMuon              ,               'mupGlobalMuon/D')
ntuple.Branch('mupTrackerMuon',            mupTrackerMuon             ,               'mupTrackerMuon/D')
# ntuple.Branch('mupStandAloneMuon',         mupStandAloneMuon          ,               'mupStandAloneMuon/D')
# ntuple.Branch('mupTMOneStationTight',      mupTMOneStationTight       ,               'mupTMOneStationTight/D')
# ntuple.Branch('mupTMOneStationLoose',      mupTMOneStationLoose       ,               'mupTMOneStationLoose/D')
# 
# ntuple.Branch('kstTrkmGlobalMuon',            kstTrkmGlobalMuon             ,               'kstTrkmGlobalMuon/D')
# ntuple.Branch('kstTrkmTrackerMuon',           kstTrkmTrackerMuon            ,               'kstTrkmTrackerMuon/D')
# ntuple.Branch('kstTrkmStandAloneMuon',        kstTrkmStandAloneMuon         ,               'kstTrkmStandAloneMuon/D')
# ntuple.Branch('kstTrkmGlobalMuonPromptTight', kstTrkmGlobalMuonPromptTight  ,               'kstTrkmGlobalMuonPromptTight/D')
# ntuple.Branch('kstTrkmTMOneStationTight',     kstTrkmTMOneStationTight      ,               'kstTrkmTMOneStationTight/D')
# ntuple.Branch('kstTrkmTMOneStationLoose',     kstTrkmTMOneStationLoose      ,               'kstTrkmTMOneStationLoose/D')
# ntuple.Branch('kstTrkmTrackerMuonArbitrated', kstTrkmTrackerMuonArbitrated  ,               'kstTrkmTrackerMuonArbitrated/D')
# ntuple.Branch('kstTrkpGlobalMuon',            kstTrkpGlobalMuon             ,               'kstTrkpGlobalMuon/D')
# ntuple.Branch('kstTrkpTrackerMuon',           kstTrkpTrackerMuon            ,               'kstTrkpTrackerMuon/D')
# ntuple.Branch('kstTrkpStandAloneMuon',        kstTrkpStandAloneMuon         ,               'kstTrkpStandAloneMuon/D')
# ntuple.Branch('kstTrkpGlobalMuonPromptTight', kstTrkpGlobalMuonPromptTight  ,               'kstTrkpGlobalMuonPromptTight/D')
# ntuple.Branch('kstTrkpTMOneStationTight',     kstTrkpTMOneStationTight      ,               'kstTrkpTMOneStationTight/D')
# ntuple.Branch('kstTrkpTMOneStationLoose',     kstTrkpTMOneStationLoose      ,               'kstTrkpTMOneStationLoose/D')
# ntuple.Branch('kstTrkpTrackerMuonArbitrated', kstTrkpTrackerMuonArbitrated  ,               'kstTrkpTrackerMuonArbitrated/D')
# 
# 
ntuple.Branch('bPt',        bPt          ,               'bPt/D')
ntuple.Branch('kstPt',      kstPt        ,               'kstPt/D')
ntuple.Branch('mumuPt',     mumuPt       ,               'mumuPt/D')
ntuple.Branch('mumPt',      mumPt        ,               'mumPt/D')
ntuple.Branch('mupPt',      mupPt        ,               'mupPt/D')
ntuple.Branch('mupCharge',  mupCharge    ,               'mupCharge/D')
ntuple.Branch('kstTrkmPt',  kstTrkmPt    ,               'kstTrkmPt/D')
ntuple.Branch('kstTrkpPt',  kstTrkpPt    ,               'kstTrkpPt/D')
ntuple.Branch('bPhi',       bPhi         ,               'bPhi/D')
ntuple.Branch('kstPhi',     kstPhi       ,               'kstPhi/D')
ntuple.Branch('mumuPhi',    mumuPhi      ,               'mumuPhi/D')
ntuple.Branch('mumPhi',     mumPhi       ,               'mumPhi/D')
ntuple.Branch('mupPhi',     mupPhi       ,               'mupPhi/D')
ntuple.Branch('kstTrkmPhi', kstTrkmPhi   ,               'kstTrkmPhi/D')
ntuple.Branch('kstTrkpPhi', kstTrkpPhi   ,               'kstTrkpPhi/D')
ntuple.Branch('bEta',       bEta         ,               'bEta/D')
ntuple.Branch('kstEta',     kstEta       ,               'kstEta/D')
ntuple.Branch('mumuEta',    mumuEta      ,               'mumuEta/D')
ntuple.Branch('mumEta',     mumEta       ,               'mumEta/D')
ntuple.Branch('mupEta',     mupEta       ,               'mupEta/D')
ntuple.Branch('kstTrkmEta', kstTrkmEta   ,               'kstTrkmEta/D')
ntuple.Branch('kstTrkpEta', kstTrkpEta   ,               'kstTrkpEta/D')
# 
# ntuple.Branch('mmk1',        mmk1          ,               'mmk1/D')
# ntuple.Branch('mmk2',        mmk2          ,               'mmk2/D')
ntuple.Branch('dR_mum_trkm', dR_mum_trkm   ,               'dR_mum_trkm/D')
ntuple.Branch('dR_mup_trkp', dR_mup_trkp   ,               'dR_mup_trkp/D')
# 
# 
# 
# ntuple.Branch('mumIsoPt_dr03'    , mumIsoPt_dr03    , 'mumIsoPt_dr03/D'   )   
# ntuple.Branch('mupIsoPt_dr03'    , mupIsoPt_dr03    , 'mupIsoPt_dr03/D'   )   
# ntuple.Branch('kstTrkmIsoPt_dr03', kstTrkmIsoPt_dr03, 'kstTrkmIsoPt_dr03/D'   )   
# ntuple.Branch('kstTrkpIsoPt_dr03', kstTrkpIsoPt_dr03, 'kstTrkpIsoPt_dr03/D'   )   
# 
# ntuple.Branch('mumIsoPt_dr04'    , mumIsoPt_dr04    , 'mumIsoPt_dr04/D'   )   
# ntuple.Branch('mupIsoPt_dr04'    , mupIsoPt_dr04    , 'mupIsoPt_dr04/D'   )   
# ntuple.Branch('kstTrkmIsoPt_dr04', kstTrkmIsoPt_dr04, 'kstTrkmIsoPt_dr04/D'   )   
# ntuple.Branch('kstTrkpIsoPt_dr04', kstTrkpIsoPt_dr04, 'kstTrkpIsoPt_dr04/D'   )   
# 



numEvents = tree_lmnr.GetEntries()
print ('total number of events in tree:', numEvents)

progressbarWidth = 40
sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
sys.stdout.flush()                          # this forces to print the stdout buffer
sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['

for i, ev in enumerate(tree_lmnr):
#     if i > 10:  continue

    if i%int(numEvents/(progressbarWidth-1))==0:
        sys.stdout.write('+')
        sys.stdout.flush()

    for var in evbr:
        var[0] = 0
    
    if not (ev.nBToTrkTrkMuMu > 0): 
      continue
      
    ## trigger requirement
    if not (ev.HLT_DoubleMu4_3_LowMass > 0 or \
            ev.HLT_DoubleMu4_3_Jpsi > 0):  
            continue  

    ## need to add here the trigger matching for the muons
    ## maybe using 
    ## TrgMatchMuon_isTriggering
    ## ??
    
    ## now loop on candidates per event
    for icand in range(ev.nBToTrkTrkMuMu):
    
        for var in newbr:
            var[0] = -99.
        for var in isobr:
            var[0] = 0.
            
        mup_index = ev.BToTrkTrkMuMu_l1_idx[icand] if ev.TrgMatchMuon_charge[ev.BToTrkTrkMuMu_l1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_l2_idx[icand] 
        mum_index = ev.BToTrkTrkMuMu_l1_idx[icand] if ev.TrgMatchMuon_charge[ev.BToTrkTrkMuMu_l1_idx[icand]] < 0 else ev.BToTrkTrkMuMu_l2_idx[icand] 
        tkp_index = ev.BToTrkTrkMuMu_trk1_idx[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_trk2_idx[icand] 
        tkm_index = ev.BToTrkTrkMuMu_trk1_idx[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] < 0 else ev.BToTrkTrkMuMu_trk2_idx[icand] 
        ditk_index = ev.BToTrkTrkMuMu_ditrack_idx[icand]

        l1_is_plus  = True if ev.TrgMatchMuon_charge[ev.BToTrkTrkMuMu_l1_idx[icand]] > 0 else False
        tk1_is_plus = True if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else False
        
        ## skim on softMuID    
        if skimSoftMu and not (ev.TrgMatchMuon_softId[mum_index] > 0 and ev.TrgMatchMuon_softId[mup_index]):
          continue

        ## per event quantities
        runN[0]                        = ev.run
#         recoVtxN[0]                    = ev.PV.npvsGood
#         trig[0]                        = paths.index(ev.TrigTable[0].split('_v')[0])
# 
#         l1_00_1p5[0]                   = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4')
#         l1_00_1p4[0]                   = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4')
#         l1_4_4[0]                      = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4')
#         l1_4p5_4p5[0]                  = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2')
# 
# #         l1_11_4[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_11_4')
# #         l1_12_5[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_12_5')
# #         l1_10_0[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_10_0_dEta_Max1p8')
# #         l1_00[0]                       = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p6_dEta_Max1p8')
# #         l1_00_OS[0]                    = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p6_dEta_Max1p8_OS')
# 
#         bsX[0]                         = ev.bsX
#         bsY[0]                         = ev.bsY
# 
        bPt[0]                         = ev.BToTrkTrkMuMu_pt[icand]
        ## why index of dimuon is not there!?
        ## for now let's build it from mu and mu2, to be changed later
        
        mumPt[0]                       = ev.BToTrkTrkMuMu_fit_l2_pt[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l1_pt[icand]
        mupPt[0]                       = ev.BToTrkTrkMuMu_fit_l1_pt[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l2_pt[icand]
        mupCharge[0]                   = ev.TrgMatchMuon_charge[mup_index]
        kstTrkmPt[0]                   = ev.BToTrkTrkMuMu_fit_trk2_pt[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk1_pt[icand]
        kstTrkpPt[0]                   = ev.BToTrkTrkMuMu_fit_trk1_pt[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk2_pt[icand]
        ## this would be the row pT
#         kstPt[0]                       = ev.DiTrack_pt[ditk_index]

        bPhi[0]                         = ev.BToTrkTrkMuMu_phi[icand]
        mumPhi[0]                       = ev.BToTrkTrkMuMu_fit_l2_phi[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l1_phi[icand]
        mupPhi[0]                       = ev.BToTrkTrkMuMu_fit_l1_phi[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l2_phi[icand]
        kstTrkmPhi[0]                   = ev.BToTrkTrkMuMu_fit_trk2_phi[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk1_phi[icand]
        kstTrkpPhi[0]                   = ev.BToTrkTrkMuMu_fit_trk1_phi[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk2_phi[icand]

        bEta[0]                         = ev.BToTrkTrkMuMu_eta[icand]
        mumEta[0]                       = ev.BToTrkTrkMuMu_fit_l2_eta[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l1_eta[icand]
        mupEta[0]                       = ev.BToTrkTrkMuMu_fit_l1_eta[icand] if l1_is_plus else ev.BToTrkTrkMuMu_fit_l2_eta[icand]
        kstTrkmEta[0]                   = ev.BToTrkTrkMuMu_fit_trk2_eta[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk1_eta[icand]
        kstTrkpEta[0]                   = ev.BToTrkTrkMuMu_fit_trk1_eta[icand] if tk1_is_plus else ev.BToTrkTrkMuMu_fit_trk2_eta[icand]

        ## need this to build the dimuon and ditrack candidate using the post fit quantities
        mum_lv_tmp  = TLorentzVector()
        mum_lv_tmp.SetPtEtaPhiM(mumPt[0], mumEta[0], mumPhi[0], muonmass_)
        mup_lv_tmp  = TLorentzVector()
        mup_lv_tmp.SetPtEtaPhiM(mupPt[0], mupEta[0], mupPhi[0], muonmass_)
        mumu_lv_tmp = mum_lv_tmp + mup_lv_tmp
        mumuPt[0] = mumu_lv_tmp.Pt()
        mumuEta[0] = mumu_lv_tmp.Eta()
        mumuPhi[0] = mumu_lv_tmp.Phi()

        ## for eta pt phi the mass assignment doesn't matter
        tkm_lv_tmp  = TLorentzVector(kstTrkmPt[0], kstTrkmEta[0], kstTrkmPhi[0], kaonmass_)
        tkp_lv_tmp  = TLorentzVector(kstTrkpPt[0], kstTrkpEta[0], kstTrkpPhi[0], pionmass_)
        kst_lv_tmp = tkm_lv_tmp + tkp_lv_tmp
        kstPt[0] = kst_lv_tmp.Pt()
        kstEta[0] = kst_lv_tmp.Eta()
        kstPhi[0] = kst_lv_tmp.Phi()
        
        

        ## per-candidate quantities
# *Br  169 :BToTrkTrkMuMu_fit_mass_Kpi :                                       *
# *         | Float_t mass of the B candidate for the leading trk->Kaon subleading trk->pion mass hypothesis*

        bMass[0]                       = ev.BToTrkTrkMuMu_fit_mass_Kpi[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_mass_piK[icand]
        bMassE[0]                      = ev.BToTrkTrkMuMu_fit_massErr_Kpi[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_massErr_piK[icand]
        bBarMass[0]                    = ev.BToTrkTrkMuMu_fit_mass_piK[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_mass_Kpi[icand]
        bBarMassE[0]                   = ev.BToTrkTrkMuMu_fit_massErr_piK[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_massErr_Kpi[icand]
        bVtxCL[0]                      = ev.BToTrkTrkMuMu_svprob[icand]  
        bVtxX[0]                       = ev.BToTrkTrkMuMu_vtx_x[icand]
        bVtxY[0]                       = ev.BToTrkTrkMuMu_vtx_y[icand]
        bVtxZ[0]                       = ev.BToTrkTrkMuMu_vtx_z[icand]
#         
        bCosAlphaBS[0]                 = ev.BToTrkTrkMuMu_fit_cos2D[icand]
#         bCosAlphaBSE[0]                = ev.bCosAlphaBSE[icand]  ## not available, not used by us later
        bLBS[0]                        = ev.BToTrkTrkMuMu_l_xy[icand]
        bLBSE[0]                       = ev.BToTrkTrkMuMu_l_xy_unc[icand]
#         bDCABS[0]                      = ev.bDCABS[icand] ### this is missing
#         bDCABSE[0]                     = ev.bDCABSE[icand] ### this is missing

        kstMass[0]                     = ev.BToTrkTrkMuMu_fit_ditrack_mass_Kpi[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_ditrack_mass_piK[icand]
#         kstMassE[0]                    = ev.DiTrack_fit_massErr_Kpi[ditk_index] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.DiTrack_fit_massErr_piK[icand]  ## missing! not used
        kstBarMass[0]                  = ev.BToTrkTrkMuMu_fit_ditrack_mass_piK[icand] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.BToTrkTrkMuMu_fit_ditrack_mass_Kpi[icand] 
#         kstBarMassE[0]                 = ev.DiTrack_fit_massErr_piK[ditk_index] if ev.Track_charge[ev.BToTrkTrkMuMu_trk1_idx[icand]] > 0 else ev.DiTrack_fit_massErr_Kpi[icand] ## missing! not used
        kstVtxCL[0]                    = ev.DiTrack_sv_prob[ditk_index]
        kstVtxX[0]                     = ev.DiTrack_vtx_x[ditk_index]
        kstVtxY[0]                     = ev.DiTrack_vtx_y[ditk_index]
        kstVtxZ[0]                     = ev.DiTrack_vtx_z[ditk_index]
        kkMass[0]                      = ev.DiTrack_fit_mass_KK[ditk_index]

        mumuMass[0]                    = mumu_lv_tmp.M()
#         mumuMassE[0]                   = ev.mumuMassE[icand]
#         mumuVtxCL[0]                   = ev.mumuVtxCL[icand]
#         mumuVtxX[0]                    = ev.mumuVtxX[icand]
#         mumuVtxY[0]                    = ev.mumuVtxY[icand]
#         mumuVtxZ[0]                    = ev.mumuVtxZ[icand]
#         mumuCosAlphaBS[0]              = ev.mumuCosAlphaBS[icand]
#         mumuCosAlphaBSE[0]             = ev.mumuCosAlphaBSE[icand]
#         mumuLBS[0]                     = ev.mumuLBS[icand]
#         mumuLBSE[0]                    = ev.mumuLBSE[icand]
#         mumuDCA[0]                     = ev.mumuDCA[icand]
#         
        ### mum_index
        mumSoft[0]                = ev.TrgMatchMuon_softId[mum_index]
        mumSoftMVA[0]             = ev.TrgMatchMuon_softMvaId[mum_index]
        mumSoftMVARun3Score[0]    = ev.TrgMatchMuon_softMvaRun3[mum_index]

        ### mup_index
        mupSoft[0]                = ev.TrgMatchMuon_softId[mup_index]      
        mupSoftMVA[0]             = ev.TrgMatchMuon_softMvaId[mup_index]   
        mupSoftMVARun3Score[0]    = ev.TrgMatchMuon_softMvaRun3[mup_index] 


#         mumHighPurity[0]               = ROOT.FindValueFromVectorOfBool(ev.mumHighPurity, icand) 
#         mumCL[0]                       = ev.mumCL[icand]
#         mumNormChi2[0]                 = ev.mumNormChi2[icand]
#         mumMinIP2D[0]                  = ev.mumMinIP2D[icand]
#         mumMinIP2DE[0]                 = ev.mumMinIP2DE[icand]
#         mumNPixLayers[0]               = ev.mumNPixLayers[icand]
#         mumNTrkLayers[0]               = ev.mumNTrkLayers[icand]
#         mumNMuonHits[0]                = ev.mumNMuonHits[icand]
#         mumNMatchStation[0]            = ev.mumNMatchStation[icand]
#         
        mumGlobalMuon            [0] = ev.TrgMatchMuon_isGlobal[mum_index]
        mumTrackerMuon           [0] = ev.TrgMatchMuon_isTracker[mum_index]
        mupGlobalMuon            [0] = ev.TrgMatchMuon_isGlobal[mup_index]
        mupTrackerMuon           [0] = ev.TrgMatchMuon_isTracker[mup_index]
#         mumStandAloneMuon        [0] = mumCategoryDict['StandAloneMuon']
#         mumTMOneStationTight     [0] = mumCategoryDict['TMOneStationTight']
#         mumTMOneStationLoose     [0] = mumCategoryDict['TMOneStationLoose']        
# 
#         mupHighPurity[0]               = ROOT.FindValueFromVectorOfBool(ev.mupHighPurity, icand)
#         mupCL[0]                       = ev.mupCL[icand]
#         mupNormChi2[0]                 = ev.mupNormChi2[icand]
#         mupdxyBS[0]                    = ev.mupdxyBS[icand]
#         mupdzBS[0]                     = ev.mupdzBS[icand]
#         mupMinIP2D[0]                  = ev.mupMinIP2D[icand]
#         mupMinIP2DE[0]                 = ev.mupMinIP2DE[icand]
#         mupNPixLayers[0]               = ev.mupNPixLayers[icand]
#         mupNTrkLayers[0]               = ev.mupNTrkLayers[icand]
#         mupNMuonHits[0]                = ev.mupNMuonHits[icand]
#         mupNMatchStation[0]            = ev.mupNMatchStation[icand]
#         
#         mupStandAloneMuon        [0] = mupCategoryDict['StandAloneMuon']
#         mupTMOneStationTight     [0] = mupCategoryDict['TMOneStationTight']
#         mupTMOneStationLoose     [0] = mupCategoryDict['TMOneStationLoose']
#         
#         if isNot2016: 
#             mumdxyBS[0]         = ev.mumdxyBS[icand]
#             mumdzBS[0]          = ev.mumdzBS[icand]
#             mupdxyBS[0]         = ev.mupdxyBS[icand]
#             mupdzBS[0]          = ev.mupdzBS[icand]
# 
#         
#         kstTrkmHighPurity[0]           = ROOT.FindValueFromVectorOfBool(ev.kstTrkmHighPurity, icand)
#         kstTrkmCL[0]                   = ev.kstTrkmCL[icand]
#         kstTrkmNormChi2[0]             = ev.kstTrkmNormChi2[icand]
#         kstTrkmPz[0]                   = ev.kstTrkmPz[icand]
#         kstTrkmDCABS[0]                = ev.kstTrkmDCABS[icand]     ## theDCAXBS.perigeeParameters().transverseImpactParameter();
#         kstTrkmDCABSE[0]               = ev.kstTrkmDCABSE[icand]
        kstTrkmDCABSSign[0]              = ev.Track_DCASig[tkm_index]
#         kstTrkmFracHits[0]             = ev.kstTrkmFracHits[icand]
        kstTrkmNPixHits[0]               = ev.Track_nValidPixelHits[tkm_index]
#         kstTrkmNPixLayers[0]           = ev.kstTrkmNPixLayers[icand]
#         kstTrkmNTrkHits[0]             = ev.kstTrkmNTrkHits[icand]
#         kstTrkmNTrkLayers[0]           = ev.kstTrkmNTrkLayers[icand]
#         kstTrkmMinIP2D[0]              = ev.kstTrkmMinIP2D[icand]
#         kstTrkmMinIP2DE[0]             = ev.kstTrkmMinIP2DE[icand]
# 
#         trkmCategoryDict = muonCategory(ev.kstTrkmMuMatch[icand])
#         kstTrkmGlobalMuon            [0] = trkmCategoryDict['GlobalMuon']
#         kstTrkmTrackerMuon           [0] = trkmCategoryDict['TrackerMuon']
#         kstTrkmStandAloneMuon        [0] = trkmCategoryDict['StandAloneMuon']
#         kstTrkmTrackerMuonArbitrated [0] = trkmCategoryDict['TrackerMuonArbitrated']
#         kstTrkmTMOneStationTight     [0] = trkmCategoryDict['TMOneStationTight']
#         kstTrkmTMOneStationLoose     [0] = trkmCategoryDict['TMOneStationLoose']
# 
#         kstTrkpHighPurity[0]           = ROOT.FindValueFromVectorOfBool(ev.kstTrkpHighPurity, icand)
#         kstTrkpCL[0]                   = ev.kstTrkpCL[icand]
#         kstTrkpNormChi2[0]             = ev.kstTrkpNormChi2[icand]
#         kstTrkpPz[0]                   = ev.kstTrkpPz[icand]
        kstTrkpDCABSSign[0]              = ev.Track_DCASig[tkp_index]
#         kstTrkpDCABSE[0]               = ev.kstTrkpDCABSE[icand]
#         kstTrkpFracHits[0]             = ev.kstTrkpFracHits[icand]
        kstTrkpNPixHits[0]               = ev.Track_nValidPixelHits[tkp_index]
#         kstTrkpNPixLayers[0]           = ev.kstTrkpNPixLayers[icand]
#         kstTrkpNTrkHits[0]             = ev.kstTrkpNTrkHits[icand]
#         kstTrkpNTrkLayers[0]           = ev.kstTrkpNTrkLayers[icand]
#         kstTrkpMinIP2D[0]              = ev.kstTrkpMinIP2D[icand]
#         kstTrkpMinIP2DE[0]             = ev.kstTrkpMinIP2DE[icand]
# 
#         trkpCategoryDict = muonCategory(ev.kstTrkpMuMatch[icand])
#         kstTrkpGlobalMuon            [0] = trkpCategoryDict['GlobalMuon']
#         kstTrkpTrackerMuon           [0] = trkpCategoryDict['TrackerMuon']
#         kstTrkpStandAloneMuon        [0] = trkpCategoryDict['StandAloneMuon']
#         kstTrkpTrackerMuonArbitrated [0] = trkpCategoryDict['TrackerMuonArbitrated']
#         kstTrkpTMOneStationTight     [0] = trkpCategoryDict['TMOneStationTight']
#         kstTrkpTMOneStationLoose     [0] = trkpCategoryDict['TMOneStationLoose']
# 
        tagB0[0]                      = FlavorTagger(kstMass[0], kstBarMass[0])  ## 1 if B0, 0 if B0bar
# 
#         bMinusVtxCL[0]               = ev.bMinusVtxCL[icand]      
#         bMinusCosAlphaBS[0]          = ev.bMinusCosAlphaBS[icand]
#         bPlusVtxCL[0]                = ev.bPlusVtxCL[icand]    
#         bPlusCosAlphaBS[0]           = ev.bPlusCosAlphaBS[icand]  
# 
# 
        cos_theta_l[0], cos_theta_k[0], phi_kst_mumu[0] = addVarsFromPtEtaPhi(
        	tagB0[0],
        	mumuPt[0],    mumuEta[0],    mumuPhi[0], mumuMass[0],
        	mumPt[0],     mumEta[0],     mumPhi[0],  
        	mupPt[0],     mupEta[0],     mupPhi[0],  
        	kstTrkmPt[0], kstTrkmEta[0], kstTrkmPhi[0],
        	kstTrkpPt[0], kstTrkpEta[0], kstTrkpPhi[0],
        	kstPt[0],     kstEta[0],     kstPhi[0],   kstMass[0], kstBarMass[0],
        	bPt[0],       bEta[0],       bPhi[0],     bMass[0],   bBarMass[0]
        )
# 
#         mmk1[0], mmk2[0] = addMMKVars(
#                           mumPt[0],      mumEta[0],      mumPhi[0],  
#                           mupPt[0],      mupEta[0],      mupPhi[0],  
#                           kstTrkpPt[0],  kstTrkpEta[0],  kstTrkpPhi[0],
#                           kstTrkmPt[0],  kstTrkmEta[0],  kstTrkmPhi[0]
#                           )
        dR_mum_trkm[0] = addDR(mumEta[0], mumPhi[0], kstTrkmEta[0], kstTrkmPhi[0])
        dR_mup_trkp[0] = addDR(mupEta[0], mupPhi[0], kstTrkpEta[0], kstTrkpPhi[0])

        if (dR_mum_trkm[0] < 1.E-4 or dR_mup_trkp[0] < 1.E-4):  continue;

#         ###########  mu - isolation ####################
        ## this is different by default
        ## I'd suggest that for validation we remove this selection and then when re-optimising the analysis we can use the isolation that is available
#         val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
#         val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
#         
#         if len( ev.mumIsoPt[icand] ) > 0:
#             for j,isoi in enumerate(ev.mumIsoPt[icand]):
#                 if ev.mumIsodR[icand][j] < 0.3:
#                     val_isoPt_dr03 += isoi
#                 if ev.mumIsodR[icand][j] < 0.4:
#                     val_isoPt_dr04 += isoi
# 
#             mumIsoPt_dr03[0] = val_isoPt_dr03
#             mumIsoPt_dr04[0] = val_isoPt_dr04
# 
#         ###########  mu + isolation ####################
#         val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
#         val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
#         
#         if len( ev.mupIsoPt[icand] ) > 0:
#             for j,isoi in enumerate(ev.mupIsoPt[icand]):
#                 if ev.mupIsodR[icand][j] < 0.3:
#                     val_isoPt_dr03 += isoi
#                 if ev.mupIsodR[icand][j] < 0.4:
#                     val_isoPt_dr04 += isoi
# 
#             mupIsoPt_dr03[0] = val_isoPt_dr03
#             mupIsoPt_dr04[0] = val_isoPt_dr04
# 
#  
#         ###########  trk - isolation ####################
#         val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
#         val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
#         
#         if len( ev.kstTrkmIsoPt[icand] ) > 0:
#             for j,isoi in enumerate(ev.kstTrkmIsoPt[icand]):
#                 if ev.kstTrkmIsodR[icand][j] < 0.3:
#                     val_isoPt_dr03 += isoi
#                 if ev.kstTrkmIsodR[icand][j] < 0.4:
#                     val_isoPt_dr04 += isoi
# 
#             kstTrkmIsoPt_dr03[0] = val_isoPt_dr03
#             kstTrkmIsoPt_dr04[0] = val_isoPt_dr04
#         ###########  trk + isolation ####################
#         val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
#         val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
# 
#         if len( ev.kstTrkpIsoPt[icand] ) > 0:
#             for j,isoi in enumerate(ev.kstTrkpIsoPt[icand]):
#                 if ev.kstTrkpIsodR[icand][j] < 0.3:
#                     val_isoPt_dr03 += isoi
#                 if ev.kstTrkpIsodR[icand][j] < 0.4:
#                     val_isoPt_dr04 += isoi
# 
#             kstTrkpIsoPt_dr03[0] = val_isoPt_dr03
#             kstTrkpIsoPt_dr04[0] = val_isoPt_dr04
# 
        ntuple. Fill()


sys.stdout.write('\n')

file_out.cd()
ntuple.Write()
file_out.Close()
    
