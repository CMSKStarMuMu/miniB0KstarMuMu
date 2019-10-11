from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-n"  , "--number"     , dest = "number"     ,  help = "number of file"              , default = '1')
parser.add_argument("-f"  , "--folder"     , dest = "folder"     ,  help = "name of output folder"       , default = 'default_folder')
# parser.add_argument("-d"  , "--dir"        , dest = "indir"      ,  help = "name of input folder (0000)" , default = '0000')
parser.add_argument("-e"  , "--era"        , dest = "theera"     ,  help = "2018B,2018B..."              , default = '2016B')
parser.add_argument("-c"  , "--chan"       , dest = "thechan"    ,  help = "LMNR, PSI"                   , default = 'LMNR' )

options = parser.parse_args()
if not options.number:   
  parser.error('Number filename not given')
  
import ROOT
ROOT.gSystem.Load('libminiB0KstarMuMuminiKstarMuMu')
from ROOT               import miniHLTObj
from DataFormats.FWLite import Events
from collections        import OrderedDict
from array              import array
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

import sys, math, itertools
sys.path.append('/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu')
from utils.angular_vars import *
from utils.utils        import *
from samples            import *

import numpy as np

###############################################
## Parameters for the flattening:

type       = options.thechan 
# eras       = []    
# eras.append(options.theera)

skim       = True
skimSoftMu = True

if type == 'LMNR':
  paths = [ 
          'HLT_DoubleMu4_LowMassNonResonantTrk_Displaced',
          ]  

if type == 'PSI':
  paths = [ 
          'HLT_DoubleMu4_JpsiTrk_Displaced',
          'HLT_DoubleMu4_PsiPrimeTrk_Displaced',
          ]  
          
###############################################
tree_lmnr = ROOT.TChain('B0KstMuMu/B0KstMuMuNTuple')

# thekey = type + '_' + options.theera
thekey = options.theera

indir = '0000'
if int(options.number) >= 1000:  indir = '0001'
if int(options.number) >= 2000:  indir = '0002'
if int(options.number) >= 3000:  indir = '0003'
if int(options.number) >= 4000:  indir = '0004'

tree_lmnr.Add('%s/%s/B0ToKstMuMu_miniaod_%s.root'%(era_dict[thekey][0],indir,options.number))

file_out = ROOT.TFile( 'ntuples/%s/ntuple_flat_%s_PostRefitMomenta_%s%s_Feb5_%s.root'%(options.folder,type,options.theera,skimSoftMu*'_skimSoftMu',options.number), 'recreate')
ntuple   = ROOT.TTree( 'ntuple', 'ntuple' )


print '@@@@@@@@@@@     flattening ntuple     @@@@@@@@@@@'
print 'input file: ', tree_lmnr.GetFile()

try:    tree_lmnr.GetEntries()
except:
  print 'ntuple tree not found'
  exit()


newbr = []
isobr = []
evbr  = []

'''
for i in branches:
    print "{NAME}\t = array('d', [-99.]);  newbr.append({NAME})".format(NAME=i).expandtabs(40)
'''
eventN                           = array('L', [ 99 ]);  
runN                             = array('d', [-99.]);  evbr .append(runN)
recoVtxN                         = array('d', [-99.]);  evbr .append(recoVtxN)
bsX                              = array('d', [-99.]);  evbr .append(bsX)
bsY                              = array('d', [-99.]);  evbr .append(bsY)

trig                             = array('d', [-99.]);  newbr.append(trig)
l1_00_1p5                        = array('d', [-99.]);  newbr.append(l1_00_1p5 )
l1_00_1p4                        = array('d', [-99.]);  newbr.append(l1_00_1p4 )
l1_4_4                           = array('d', [-99.]);  newbr.append(l1_4_4    )
l1_4p5_4p5                       = array('d', [-99.]);  newbr.append(l1_4p5_4p5)
# l1_11_4                          = array('d', [-99.]);  newbr.append(l1_11_4 )
# l1_12_5                          = array('d', [-99.]);  newbr.append(l1_12_5 )
# l1_10_0                          = array('d', [-99.]);  newbr.append(l1_10_0 )
# l1_00                            = array('d', [-99.]);  newbr.append(l1_00   )
# l1_00_OS                         = array('d', [-99.]);  newbr.append(l1_00_OS)
bMass                            = array('d', [-99.]);  newbr.append(bMass)
bMassE                           = array('d', [-99.]);  newbr.append(bMassE)
bBarMass                         = array('d', [-99.]);  newbr.append(bBarMass)
bBarMassE                        = array('d', [-99.]);  newbr.append(bBarMassE)
bVtxCL                           = array('d', [-99.]);  newbr.append(bVtxCL)
bVtxX                            = array('d', [-99.]);  newbr.append(bVtxX)
bVtxY                            = array('d', [-99.]);  newbr.append(bVtxY)
bVtxZ                            = array('d', [-99.]);  newbr.append(bVtxZ)
bCosAlphaBS                      = array('d', [-99.]);  newbr.append(bCosAlphaBS)
bCosAlphaBSE                     = array('d', [-99.]);  newbr.append(bCosAlphaBSE)
bLBS                             = array('d', [-99.]);  newbr.append(bLBS)
bLBSE                            = array('d', [-99.]);  newbr.append(bLBSE)
bDCABS                           = array('d', [-99.]);  newbr.append(bDCABS)
bDCABSE                          = array('d', [-99.]);  newbr.append(bDCABSE)

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
mumuMassE                        = array('d', [-99.]);  newbr.append(mumuMassE)
mumuVtxCL                        = array('d', [-99.]);  newbr.append(mumuVtxCL)
mumuVtxX                         = array('d', [-99.]);  newbr.append(mumuVtxX)
mumuVtxY                         = array('d', [-99.]);  newbr.append(mumuVtxY)
mumuVtxZ                         = array('d', [-99.]);  newbr.append(mumuVtxZ)
mumuCosAlphaBS                   = array('d', [-99.]);  newbr.append(mumuCosAlphaBS)
mumuCosAlphaBSE                  = array('d', [-99.]);  newbr.append(mumuCosAlphaBSE)
mumuLBS                          = array('d', [-99.]);  newbr.append(mumuLBS)
mumuLBSE                         = array('d', [-99.]);  newbr.append(mumuLBSE)
mumuDCA                          = array('d', [-99.]);  newbr.append(mumuDCA)

mumHighPurity                    = array('d', [-99.]);  newbr.append(mumHighPurity)
mumCL                            = array('d', [-99.]);  newbr.append(mumCL)
mumNormChi2                      = array('d', [-99.]);  newbr.append(mumNormChi2)
mumDCAVtx                        = array('d', [-99.]);  newbr.append(mumDCAVtx)
mumDCAVtxE                       = array('d', [-99.]);  newbr.append(mumDCAVtxE)
mumDCABS                         = array('d', [-99.]);  newbr.append(mumDCABS)
mumDCABSE                        = array('d', [-99.]);  newbr.append(mumDCABSE)
mumdxyBS                         = array('d', [-99.]);  newbr.append(mumdxyBS)
mumdzBS                          = array('d', [-99.]);  newbr.append(mumdzBS)
mumMinIP2D                       = array('d', [-99.]);  newbr.append(mumMinIP2D)
mumMinIP2DE                      = array('d', [-99.]);  newbr.append(mumMinIP2DE)
mumNPixLayers                    = array('d', [-99.]);  newbr.append(mumNPixLayers)
mumNTrkLayers                    = array('d', [-99.]);  newbr.append(mumNTrkLayers)
mumNMuonHits                     = array('d', [-99.]);  newbr.append(mumNMuonHits)
mumNMatchStation                 = array('d', [-99.]);  newbr.append(mumNMatchStation)
mupHighPurity                    = array('d', [-99.]);  newbr.append(mupHighPurity)
mupCL                            = array('d', [-99.]);  newbr.append(mupCL)
mupNormChi2                      = array('d', [-99.]);  newbr.append(mupNormChi2)
mupDCAVtx                        = array('d', [-99.]);  newbr.append(mupDCAVtx)
mupDCAVtxE                       = array('d', [-99.]);  newbr.append(mupDCAVtxE)
mupDCABS                         = array('d', [-99.]);  newbr.append(mupDCABS)
mupDCABSE                        = array('d', [-99.]);  newbr.append(mupDCABSE)
mupdxyBS                         = array('d', [-99.]);  newbr.append(mupdxyBS)
mupdzBS                          = array('d', [-99.]);  newbr.append(mupdzBS)
mupMinIP2D                       = array('d', [-99.]);  newbr.append(mupMinIP2D)
mupMinIP2DE                      = array('d', [-99.]);  newbr.append(mupMinIP2DE)
mupNPixLayers                    = array('d', [-99.]);  newbr.append(mupNPixLayers)
mupNTrkLayers                    = array('d', [-99.]);  newbr.append(mupNTrkLayers)
mupNMuonHits                     = array('d', [-99.]);  newbr.append(mupNMuonHits)
mupNMatchStation                 = array('d', [-99.]);  newbr.append(mupNMatchStation)
kstTrkmHighPurity                = array('d', [-99.]);  newbr.append(kstTrkmHighPurity)
kstTrkmCL                        = array('d', [-99.]);  newbr.append(kstTrkmCL)
kstTrkmNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkmNormChi2)
kstTrkmPz                        = array('d', [-99.]);  newbr.append(kstTrkmPz)
kstTrkmDCABS                     = array('d', [-99.]);  newbr.append(kstTrkmDCABS)
kstTrkmDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkmDCABSE)
kstTrkmFracHits                  = array('d', [-99.]);  newbr.append(kstTrkmFracHits)
kstTrkmMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkmMinIP2D)
kstTrkmMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkmMinIP2DE)
kstTrkmNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkmNPixHits)
kstTrkmNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkmNPixLayers)
kstTrkmNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkmNTrkHits)
kstTrkmNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkmNTrkLayers)
kstTrkmMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkmMuMatch)
kstTrkpHighPurity                = array('d', [-99.]);  newbr.append(kstTrkpHighPurity)
kstTrkpCL                        = array('d', [-99.]);  newbr.append(kstTrkpCL)
kstTrkpNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkpNormChi2)
kstTrkpPz                        = array('d', [-99.]);  newbr.append(kstTrkpPz)
kstTrkpDCABS                     = array('d', [-99.]);  newbr.append(kstTrkpDCABS)
kstTrkpDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkpDCABSE)
kstTrkpFracHits                  = array('d', [-99.]);  newbr.append(kstTrkpFracHits)
kstTrkpMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkpMinIP2D)
kstTrkpMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkpMinIP2DE)
kstTrkpNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkpNPixHits)
kstTrkpNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkpNPixLayers)
kstTrkpNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkpNTrkHits)
kstTrkpNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkpNTrkLayers)
kstTrkpMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkpMuMatch)
tagB0                            = array('d', [-99.]);  newbr.append(tagB0)

bMinusVtxCL                      = array('d', [-99.]);  newbr.append(bMinusVtxCL      )
bMinusCosAlphaBS                 = array('d', [-99.]);  newbr.append(bMinusCosAlphaBS )
bPlusVtxCL                       = array('d', [-99.]);  newbr.append(bPlusVtxCL       )
bPlusCosAlphaBS                  = array('d', [-99.]);  newbr.append(bPlusCosAlphaBS  )

mumGlobalMuon                    = array('d',[-99.]); newbr.append(mumGlobalMuon             )   
mumTrackerMuon                   = array('d',[-99.]); newbr.append(mumTrackerMuon            )   
mumStandAloneMuon                = array('d',[-99.]); newbr.append(mumStandAloneMuon         )   
mumTMOneStationTight             = array('d',[-99.]); newbr.append(mumTMOneStationTight      )   
mumTMOneStationLoose             = array('d',[-99.]); newbr.append(mumTMOneStationLoose      )   
mupGlobalMuon                    = array('d',[-99.]); newbr.append(mupGlobalMuon             )   
mupTrackerMuon                   = array('d',[-99.]); newbr.append(mupTrackerMuon            )   
mupStandAloneMuon                = array('d',[-99.]); newbr.append(mupStandAloneMuon         )   
mupTMOneStationTight             = array('d',[-99.]); newbr.append(mupTMOneStationTight      )   
mupTMOneStationLoose             = array('d',[-99.]); newbr.append(mupTMOneStationLoose      )   
kstTrkmGlobalMuon                = array('d',[-99.]); newbr.append(kstTrkmGlobalMuon             )   
kstTrkmTrackerMuon               = array('d',[-99.]); newbr.append(kstTrkmTrackerMuon            )   
kstTrkmStandAloneMuon            = array('d',[-99.]); newbr.append(kstTrkmStandAloneMuon         )   
kstTrkmGlobalMuonPromptTight     = array('d',[-99.]); newbr.append(kstTrkmGlobalMuonPromptTight  )   
kstTrkmTMOneStationTight         = array('d',[-99.]); newbr.append(kstTrkmTMOneStationTight      )   
kstTrkmTMOneStationLoose         = array('d',[-99.]); newbr.append(kstTrkmTMOneStationLoose      )   
kstTrkmTrackerMuonArbitrated     = array('d',[-99.]); newbr.append(kstTrkmTrackerMuonArbitrated  )   
kstTrkpGlobalMuon                = array('d',[-99.]); newbr.append(kstTrkpGlobalMuon             )   
kstTrkpTrackerMuon               = array('d',[-99.]); newbr.append(kstTrkpTrackerMuon            )   
kstTrkpStandAloneMuon            = array('d',[-99.]); newbr.append(kstTrkpStandAloneMuon         )   
kstTrkpGlobalMuonPromptTight     = array('d',[-99.]); newbr.append(kstTrkpGlobalMuonPromptTight  )   
kstTrkpTMOneStationTight         = array('d',[-99.]); newbr.append(kstTrkpTMOneStationTight      )   
kstTrkpTMOneStationLoose         = array('d',[-99.]); newbr.append(kstTrkpTMOneStationLoose      )   
kstTrkpTrackerMuonArbitrated     = array('d',[-99.]); newbr.append(kstTrkpTrackerMuonArbitrated  )   


bPt                              = array('d',[-99.]); newbr.append( bPt          )   
kstPt                            = array('d',[-99.]); newbr.append( kstPt        )   
mumuPt                           = array('d',[-99.]); newbr.append( mumuPt       )   
mumPt                            = array('d',[-99.]); newbr.append( mumPt        )   
mupPt                            = array('d',[-99.]); newbr.append( mupPt        )   
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

mmk1                              = array('d',[-99.]); newbr.append( mmk1          )   
mmk2                              = array('d',[-99.]); newbr.append( mmk2          )   
dR_mum_trkm                       = array('d',[-99.]); newbr.append( dR_mum_trkm   )          
dR_mup_trkp                       = array('d',[-99.]); newbr.append( dR_mup_trkp   )          


charge_trig_matched              = array('d',[0.])   ; newbr.append( charge_trig_matched   )   
cos_theta_l                      = array('d', [-99.]); newbr.append(cos_theta_l)
cos_theta_k                      = array('d', [-99.]); newbr.append(cos_theta_k)
phi_kst_mumu                     = array('d', [-99.]); newbr.append(phi_kst_mumu)

mumIsoPt_dr03                    = array('d',[0.]); isobr.append( mumIsoPt_dr03   )   
mupIsoPt_dr03                    = array('d',[0.]); isobr.append( mupIsoPt_dr03   )   
kstTrkmIsoPt_dr03                = array('d',[0.]); isobr.append( kstTrkmIsoPt_dr03   )   
kstTrkpIsoPt_dr03                = array('d',[0.]); isobr.append( kstTrkpIsoPt_dr03   )   

mumIsoPt_dr04                    = array('d',[0.]); isobr.append( mumIsoPt_dr04   )   
mupIsoPt_dr04                    = array('d',[0.]); isobr.append( mupIsoPt_dr04   )   
kstTrkmIsoPt_dr04                = array('d',[0.]); isobr.append( kstTrkmIsoPt_dr04   )   
kstTrkpIsoPt_dr04                = array('d',[0.]); isobr.append( kstTrkpIsoPt_dr04   )   

mumIsoPt_dr05                    = array('d',[0.]); isobr.append( mumIsoPt_dr05   )   
mupIsoPt_dr05                    = array('d',[0.]); isobr.append( mupIsoPt_dr05   )   
kstTrkmIsoPt_dr05                = array('d',[0.]); isobr.append( kstTrkmIsoPt_dr05   )   
kstTrkpIsoPt_dr05                = array('d',[0.]); isobr.append( kstTrkpIsoPt_dr05   )   




if not skim:
  mumNPixHits                      = array('d', [-99.]);  newbr.append(mumNPixHits)
  mumNTrkHits                      = array('d', [-99.]);  newbr.append(mumNTrkHits)
  mupNPixHits                      = array('d', [-99.]);  newbr.append(mupNPixHits)
  mupNTrkHits                      = array('d', [-99.]);  newbr.append(mupNTrkHits)
  mumKinkChi2                      = array('d', [-99.]);  newbr.append(mumKinkChi2)
  mumFracHits                      = array('d', [-99.]);  newbr.append(mumFracHits)
  mupKinkChi2                      = array('d', [-99.]);  newbr.append(mupKinkChi2)
  mupFracHits                      = array('d', [-99.]);  newbr.append(mupFracHits)
  mumMinIP                         = array('d', [-99.]);  newbr.append(mumMinIP)
  mumMinIPS                        = array('d', [-99.]);  newbr.append(mumMinIPS)
  mupMinIP                         = array('d', [-99.]);  newbr.append(mupMinIP)
  mupMinIPS                        = array('d', [-99.]);  newbr.append(mupMinIPS)
  kstTrkmDCAVtx                    = array('d', [-99.]);  newbr.append(kstTrkmDCAVtx)
  kstTrkmDCAVtxE                   = array('d', [-99.]);  newbr.append(kstTrkmDCAVtxE)
  kstTrkpDCAVtx                    = array('d', [-99.]);  newbr.append(kstTrkpDCAVtx)
  kstTrkpDCAVtxE                   = array('d', [-99.]);  newbr.append(kstTrkpDCAVtxE)
  mumGlobalMuonPromptTight         = array('d',[-99.]); newbr.append(mumGlobalMuonPromptTight  )   
  mumTMLastStationTight            = array('d',[-99.]); newbr.append(mumTMLastStationTight     )   
  mumTMLastStationLoose            = array('d',[-99.]); newbr.append(mumTMLastStationLoose     )   
  mumTM2DCompatibilityTight        = array('d',[-99.]); newbr.append(mumTM2DCompatibilityTight )   
  mumTM2DCompatibilityLoose        = array('d',[-99.]); newbr.append(mumTM2DCompatibilityLoose )   
  mumTMLastStationAngTight         = array('d',[-99.]); newbr.append(mumTMLastStationAngTight  )   
  mumTMLastStationAngLoose         = array('d',[-99.]); newbr.append(mumTMLastStationAngLoose  )   
  mumTMOneStationAngTight          = array('d',[-99.]); newbr.append(mumTMOneStationAngTight   )   
  mumTMOneStationAngLoose          = array('d',[-99.]); newbr.append(mumTMOneStationAngLoose   )   
  mumTrackerMuonArbitrated         = array('d',[-99.]); newbr.append(mumTrackerMuonArbitrated  )   
  mupGlobalMuonPromptTight         = array('d',[-99.]); newbr.append(mupGlobalMuonPromptTight  )   
  mupTMLastStationTight            = array('d',[-99.]); newbr.append(mupTMLastStationTight     )   
  mupTMLastStationLoose            = array('d',[-99.]); newbr.append(mupTMLastStationLoose     )   
  mupTM2DCompatibilityTight        = array('d',[-99.]); newbr.append(mupTM2DCompatibilityTight )   
  mupTM2DCompatibilityLoose        = array('d',[-99.]); newbr.append(mupTM2DCompatibilityLoose )   
  mupTMLastStationAngTight         = array('d',[-99.]); newbr.append(mupTMLastStationAngTight  )   
  mupTMLastStationAngLoose         = array('d',[-99.]); newbr.append(mupTMLastStationAngLoose  )   
  mupTMOneStationAngTight          = array('d',[-99.]); newbr.append(mupTMOneStationAngTight   )   
  mupTMOneStationAngLoose          = array('d',[-99.]); newbr.append(mupTMOneStationAngLoose   )   
  mupTrackerMuonArbitrated         = array('d',[-99.]); newbr.append(mupTrackerMuonArbitrated  )   
  kstTrkmTM2DCompatibilityTight    = array('d',[-99.]); newbr.append(kstTrkmTM2DCompatibilityTight )   
  kstTrkmTM2DCompatibilityLoose    = array('d',[-99.]); newbr.append(kstTrkmTM2DCompatibilityLoose )   
  kstTrkpTM2DCompatibilityTight    = array('d',[-99.]); newbr.append(kstTrkpTM2DCompatibilityTight )   
  kstTrkpTM2DCompatibilityLoose    = array('d',[-99.]); newbr.append(kstTrkpTM2DCompatibilityLoose )   
  kstTrkmTMLastStationAngTight     = array('d',[-99.]); newbr.append(kstTrkmTMLastStationAngTight  )   
  kstTrkmTMLastStationAngLoose     = array('d',[-99.]); newbr.append(kstTrkmTMLastStationAngLoose  )   
  kstTrkmTMOneStationAngTight      = array('d',[-99.]); newbr.append(kstTrkmTMOneStationAngTight   )   
  kstTrkmTMOneStationAngLoose      = array('d',[-99.]); newbr.append(kstTrkmTMOneStationAngLoose   )   
  kstTrkpTMLastStationAngTight     = array('d',[-99.]); newbr.append(kstTrkpTMLastStationAngTight  )   
  kstTrkpTMLastStationAngLoose     = array('d',[-99.]); newbr.append(kstTrkpTMLastStationAngLoose  )   
  kstTrkpTMOneStationAngTight      = array('d',[-99.]); newbr.append(kstTrkpTMOneStationAngTight   )   
  kstTrkpTMOneStationAngLoose      = array('d',[-99.]); newbr.append(kstTrkpTMOneStationAngLoose   )   
  kstTrkmTMLastStationTight        = array('d',[-99.]); newbr.append(kstTrkmTMLastStationTight     )   
  kstTrkmTMLastStationLoose        = array('d',[-99.]); newbr.append(kstTrkmTMLastStationLoose     )   
  kstTrkpTMLastStationTight        = array('d',[-99.]); newbr.append(kstTrkpTMLastStationTight     )   
  kstTrkpTMLastStationLoose        = array('d',[-99.]); newbr.append(kstTrkpTMLastStationLoose     )   
  kstTrkmMinIP                     = array('d',[-99.]);  newbr.append(kstTrkmMinIP)
  kstTrkmMinIPS                    = array('d',[-99.]);  newbr.append(kstTrkmMinIPS)
  kstTrkpMinIP                     = array('d',[-99.]);  newbr.append(kstTrkpMinIP)
  kstTrkpMinIPS                    = array('d',[-99.]);  newbr.append(kstTrkpMinIPS)

  mumIso                           = ROOT.std.vector('double')()  
  mupIso                           = ROOT.std.vector('double')()  
  kstTrkmIso                       = ROOT.std.vector('double')()  
  kstTrkpIso                       = ROOT.std.vector('double')()  
    
    
    
'''
for i in branches:
    print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/D')".format(NAME=i).expandtabs(40)
'''
ntuple.Branch('runN',                   runN,                                   'runN/D')
ntuple.Branch('eventN',                 eventN,                                 'eventN/L')
ntuple.Branch('recoVtxN',               recoVtxN,                               'recoVtxN/D')
ntuple.Branch('bsX',                    bsX,                                    'bsX/D')
ntuple.Branch('bsY',                    bsY,                                    'bsY/D')

ntuple.Branch('trig',                   trig,                                   'trig/D')
ntuple.Branch('l1_00_1p5',              l1_00_1p5 ,                             'l1_00_1p5/D')
ntuple.Branch('l1_00_1p4',              l1_00_1p4 ,                             'l1_00_1p4/D')
ntuple.Branch('l1_4_4',                 l1_4_4 ,                                'l1_4_4/D')
ntuple.Branch('l1_4p5_4p5',             l1_4p5_4p5 ,                            'l1_4p5_4p5/D')

ntuple.Branch('bMass',                  bMass,                                  'bMass/D')
ntuple.Branch('bMassE',                 bMassE,                                 'bMassE/D')
ntuple.Branch('bBarMass',               bBarMass,                               'bBarMass/D')
ntuple.Branch('bBarMassE',              bBarMassE,                              'bBarMassE/D')
ntuple.Branch('bVtxCL',                 bVtxCL,                                 'bVtxCL/D')
ntuple.Branch('bVtxX',                  bVtxX,                                  'bVtxX/D')
ntuple.Branch('bVtxY',                  bVtxY,                                  'bVtxY/D')
ntuple.Branch('bVtxZ',                  bVtxZ,                                  'bVtxZ/D')

ntuple.Branch('bCosAlphaBS',            bCosAlphaBS,                            'bCosAlphaBS/D')
ntuple.Branch('bCosAlphaBSE',           bCosAlphaBSE,                           'bCosAlphaBSE/D')
ntuple.Branch('bLBS',                   bLBS,                                   'bLBS/D')
ntuple.Branch('bLBSE',                  bLBSE,                                  'bLBSE/D')
ntuple.Branch('bDCABS',                 bDCABS,                                 'bDCABS/D')
ntuple.Branch('bDCABSE',                bDCABSE,                                'bDCABSE/D')

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
ntuple.Branch('mumuMassE',              mumuMassE,                              'mumuMassE/D')
ntuple.Branch('mumuVtxCL',              mumuVtxCL,                              'mumuVtxCL/D')
ntuple.Branch('mumuVtxX',               mumuVtxX,                               'mumuVtxX/D')
ntuple.Branch('mumuVtxY',               mumuVtxY,                               'mumuVtxY/D')
ntuple.Branch('mumuVtxZ',               mumuVtxZ,                               'mumuVtxZ/D')
ntuple.Branch('mumuCosAlphaBS',         mumuCosAlphaBS,                         'mumuCosAlphaBS/D')
ntuple.Branch('mumuCosAlphaBSE',        mumuCosAlphaBSE,                        'mumuCosAlphaBSE/D')
ntuple.Branch('mumuLBS',                mumuLBS,                                'mumuLBS/D')
ntuple.Branch('mumuLBSE',               mumuLBSE,                               'mumuLBSE/D')
ntuple.Branch('mumuDCA',                mumuDCA,                                'mumuDCA/D')

ntuple.Branch('mumHighPurity',          mumHighPurity,                          'mumHighPurity/D')
ntuple.Branch('mumCL',                  mumCL,                                  'mumCL/D')
ntuple.Branch('mumNormChi2',            mumNormChi2,                            'mumNormChi2/D')
ntuple.Branch('mumDCAVtx',              mumDCAVtx,                              'mumDCAVtx/D')
ntuple.Branch('mumDCAVtxE',             mumDCAVtxE,                             'mumDCAVtxE/D')
ntuple.Branch('mumDCABS',               mumDCABS,                               'mumDCABS/D')
ntuple.Branch('mumDCABSE',              mumDCABSE,                              'mumDCABSE/D')
ntuple.Branch('mumdxyBS',               mumdxyBS,                               'mumdxyBS/D')
ntuple.Branch('mumdzBS',                mumdzBS,                                'mumdzBS/D')
ntuple.Branch('mumMinIP2D',             mumMinIP2D,                             'mumMinIP2D/D')
ntuple.Branch('mumMinIP2DE',            mumMinIP2DE,                            'mumMinIP2DE/D')
ntuple.Branch('mumNPixLayers',          mumNPixLayers,                          'mumNPixLayers/D')
ntuple.Branch('mumNTrkLayers',          mumNTrkLayers,                          'mumNTrkLayers/D')
ntuple.Branch('mumNMuonHits',           mumNMuonHits,                           'mumNMuonHits/D')
ntuple.Branch('mumNMatchStation',       mumNMatchStation,                       'mumNMatchStation/D')
ntuple.Branch('mupHighPurity',          mupHighPurity,                          'mupHighPurity/D')
ntuple.Branch('mupCL',                  mupCL,                                  'mupCL/D')
ntuple.Branch('mupNormChi2',            mupNormChi2,                            'mupNormChi2/D')
ntuple.Branch('mupDCAVtx',              mupDCAVtx,                              'mupDCAVtx/D')
ntuple.Branch('mupDCAVtxE',             mupDCAVtxE,                             'mupDCAVtxE/D')
ntuple.Branch('mupDCABS',               mupDCABS,                               'mupDCABS/D')
ntuple.Branch('mupDCABSE',              mupDCABSE,                              'mupDCABSE/D')
ntuple.Branch('mupdxyBS',               mupdxyBS,                               'mupdxyBS/D')
ntuple.Branch('mupdzBS',                mupdzBS,                                'mupdzBS/D')
ntuple.Branch('mupMinIP2D',             mupMinIP2D,                             'mupMinIP2D/D')
ntuple.Branch('mupMinIP2DE',            mupMinIP2DE,                            'mupMinIP2DE/D')
ntuple.Branch('mupNPixLayers',          mupNPixLayers,                          'mupNPixLayers/D')
ntuple.Branch('mupNTrkLayers',          mupNTrkLayers,                          'mupNTrkLayers/D')
ntuple.Branch('mupNMuonHits',           mupNMuonHits,                           'mupNMuonHits/D')
ntuple.Branch('mupNMatchStation',       mupNMatchStation,                       'mupNMatchStation/D')

ntuple.Branch('kstTrkmHighPurity',      kstTrkmHighPurity,                      'kstTrkmHighPurity/D')
ntuple.Branch('kstTrkmCL',              kstTrkmCL,                              'kstTrkmCL/D')
ntuple.Branch('kstTrkmNormChi2',        kstTrkmNormChi2,                        'kstTrkmNormChi2/D')
ntuple.Branch('kstTrkmPz',              kstTrkmPz,                              'kstTrkmPz/D')
ntuple.Branch('kstTrkmDCABS',           kstTrkmDCABS,                           'kstTrkmDCABS/D')
ntuple.Branch('kstTrkmDCABSE',          kstTrkmDCABSE,                          'kstTrkmDCABSE/D')
ntuple.Branch('kstTrkmFracHits',        kstTrkmFracHits,                        'kstTrkmFracHits/D')
ntuple.Branch('kstTrkmMinIP2D',         kstTrkmMinIP2D,                         'kstTrkmMinIP2D/D')
ntuple.Branch('kstTrkmMinIP2DE',        kstTrkmMinIP2DE,                        'kstTrkmMinIP2DE/D')
ntuple.Branch('kstTrkmNPixHits',        kstTrkmNPixHits,                        'kstTrkmNPixHits/D')
ntuple.Branch('kstTrkmNPixLayers',      kstTrkmNPixLayers,                      'kstTrkmNPixLayers/D')
ntuple.Branch('kstTrkmNTrkHits',        kstTrkmNTrkHits,                        'kstTrkmNTrkHits/D')
ntuple.Branch('kstTrkmNTrkLayers',      kstTrkmNTrkLayers,                      'kstTrkmNTrkLayers/D')
ntuple.Branch('kstTrkmMuMatch',         kstTrkmMuMatch,                         'kstTrkmMuMatch/D')
ntuple.Branch('kstTrkpHighPurity',      kstTrkpHighPurity,                      'kstTrkpHighPurity/D')
ntuple.Branch('kstTrkpCL',              kstTrkpCL,                              'kstTrkpCL/D')
ntuple.Branch('kstTrkpNormChi2',        kstTrkpNormChi2,                        'kstTrkpNormChi2/D')
ntuple.Branch('kstTrkpPz',              kstTrkpPz,                              'kstTrkpPz/D')
ntuple.Branch('kstTrkpDCABS',           kstTrkpDCABS,                           'kstTrkpDCABS/D')
ntuple.Branch('kstTrkpDCABSE',          kstTrkpDCABSE,                          'kstTrkpDCABSE/D')
ntuple.Branch('kstTrkpFracHits',        kstTrkpFracHits,                        'kstTrkpFracHits/D')
ntuple.Branch('kstTrkpMinIP2D',         kstTrkpMinIP2D,                         'kstTrkpMinIP2D/D')
ntuple.Branch('kstTrkpMinIP2DE',        kstTrkpMinIP2DE,                        'kstTrkpMinIP2DE/D')
ntuple.Branch('kstTrkpNPixHits',        kstTrkpNPixHits,                        'kstTrkpNPixHits/D')
ntuple.Branch('kstTrkpNPixLayers',      kstTrkpNPixLayers,                      'kstTrkpNPixLayers/D')
ntuple.Branch('kstTrkpNTrkHits',        kstTrkpNTrkHits,                        'kstTrkpNTrkHits/D')
ntuple.Branch('kstTrkpNTrkLayers',      kstTrkpNTrkLayers,                      'kstTrkpNTrkLayers/D')
ntuple.Branch('kstTrkpMuMatch',         kstTrkpMuMatch,                         'kstTrkpMuMatch/D')
ntuple.Branch('tagB0',                  tagB0,                                  'tagB0/D')
ntuple.Branch('cos_theta_l' ,           cos_theta_l ,                           'cos_theta_l/D')
ntuple.Branch('cos_theta_k' ,           cos_theta_k ,                           'cos_theta_k/D')
ntuple.Branch('phi_kst_mumu',           phi_kst_mumu,                           'phi_kst_mumu/D')

ntuple.Branch('bMinusVtxCL',            bMinusVtxCL     ,                       'bMinusVtxCL/D')
ntuple.Branch('bMinusCosAlphaBS',       bMinusCosAlphaBS,                       'bMinusCosAlphaBS/D')
ntuple.Branch('bPlusVtxCL',             bPlusVtxCL      ,                       'bPlusVtxCL/D')
ntuple.Branch('bPlusCosAlphaBS',        bPlusCosAlphaBS ,                       'bPlusCosAlphaBS/D')

ntuple.Branch('mumGlobalMuon',             mumGlobalMuon              ,               'mumGlobalMuon/D')
ntuple.Branch('mumTrackerMuon',            mumTrackerMuon             ,               'mumTrackerMuon/D')
ntuple.Branch('mumStandAloneMuon',         mumStandAloneMuon          ,               'mumStandAloneMuon/D')
ntuple.Branch('mumTMOneStationTight',      mumTMOneStationTight       ,               'mumTMOneStationTight/D')
ntuple.Branch('mumTMOneStationLoose',      mumTMOneStationLoose       ,               'mumTMOneStationLoose/D')
ntuple.Branch('mupGlobalMuon',             mupGlobalMuon              ,               'mupGlobalMuon/D')
ntuple.Branch('mupTrackerMuon',            mupTrackerMuon             ,               'mupTrackerMuon/D')
ntuple.Branch('mupStandAloneMuon',         mupStandAloneMuon          ,               'mupStandAloneMuon/D')
ntuple.Branch('mupTMOneStationTight',      mupTMOneStationTight       ,               'mupTMOneStationTight/D')
ntuple.Branch('mupTMOneStationLoose',      mupTMOneStationLoose       ,               'mupTMOneStationLoose/D')

ntuple.Branch('kstTrkmGlobalMuon',            kstTrkmGlobalMuon             ,               'kstTrkmGlobalMuon/D')
ntuple.Branch('kstTrkmTrackerMuon',           kstTrkmTrackerMuon            ,               'kstTrkmTrackerMuon/D')
ntuple.Branch('kstTrkmStandAloneMuon',        kstTrkmStandAloneMuon         ,               'kstTrkmStandAloneMuon/D')
ntuple.Branch('kstTrkmGlobalMuonPromptTight', kstTrkmGlobalMuonPromptTight  ,               'kstTrkmGlobalMuonPromptTight/D')
ntuple.Branch('kstTrkmTMOneStationTight',     kstTrkmTMOneStationTight      ,               'kstTrkmTMOneStationTight/D')
ntuple.Branch('kstTrkmTMOneStationLoose',     kstTrkmTMOneStationLoose      ,               'kstTrkmTMOneStationLoose/D')
ntuple.Branch('kstTrkmTrackerMuonArbitrated', kstTrkmTrackerMuonArbitrated  ,               'kstTrkmTrackerMuonArbitrated/D')
ntuple.Branch('kstTrkpGlobalMuon',            kstTrkpGlobalMuon             ,               'kstTrkpGlobalMuon/D')
ntuple.Branch('kstTrkpTrackerMuon',           kstTrkpTrackerMuon            ,               'kstTrkpTrackerMuon/D')
ntuple.Branch('kstTrkpStandAloneMuon',        kstTrkpStandAloneMuon         ,               'kstTrkpStandAloneMuon/D')
ntuple.Branch('kstTrkpGlobalMuonPromptTight', kstTrkpGlobalMuonPromptTight  ,               'kstTrkpGlobalMuonPromptTight/D')
ntuple.Branch('kstTrkpTMOneStationTight',     kstTrkpTMOneStationTight      ,               'kstTrkpTMOneStationTight/D')
ntuple.Branch('kstTrkpTMOneStationLoose',     kstTrkpTMOneStationLoose      ,               'kstTrkpTMOneStationLoose/D')
ntuple.Branch('kstTrkpTrackerMuonArbitrated', kstTrkpTrackerMuonArbitrated  ,               'kstTrkpTrackerMuonArbitrated/D')


ntuple.Branch('bPt',        bPt          ,               'bPt/D')
ntuple.Branch('kstPt',      kstPt        ,               'kstPt/D')
ntuple.Branch('mumuPt',     mumuPt       ,               'mumuPt/D')
ntuple.Branch('mumPt',      mumPt        ,               'mumPt/D')
ntuple.Branch('mupPt',      mupPt        ,               'mupPt/D')
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

ntuple.Branch('mmk1',        mmk1          ,               'mmk1/D')
ntuple.Branch('mmk2',        mmk2          ,               'mmk2/D')
ntuple.Branch('dR_mum_trkm', dR_mum_trkm   ,               'dR_mum_trkm/D')
ntuple.Branch('dR_mup_trkp', dR_mup_trkp   ,               'dR_mup_trkp/D')



ntuple.Branch('mumIsoPt_dr03'    , mumIsoPt_dr03    , 'mumIsoPt_dr03/D'   )   
ntuple.Branch('mupIsoPt_dr03'    , mupIsoPt_dr03    , 'mupIsoPt_dr03/D'   )   
ntuple.Branch('kstTrkmIsoPt_dr03', kstTrkmIsoPt_dr03, 'kstTrkmIsoPt_dr03/D'   )   
ntuple.Branch('kstTrkpIsoPt_dr03', kstTrkpIsoPt_dr03, 'kstTrkpIsoPt_dr03/D'   )   

ntuple.Branch('mumIsoPt_dr04'    , mumIsoPt_dr04    , 'mumIsoPt_dr04/D'   )   
ntuple.Branch('mupIsoPt_dr04'    , mupIsoPt_dr04    , 'mupIsoPt_dr04/D'   )   
ntuple.Branch('kstTrkmIsoPt_dr04', kstTrkmIsoPt_dr04, 'kstTrkmIsoPt_dr04/D'   )   
ntuple.Branch('kstTrkpIsoPt_dr04', kstTrkpIsoPt_dr04, 'kstTrkpIsoPt_dr04/D'   )   

ntuple.Branch('mumIsoPt_dr05'    , mumIsoPt_dr05    , 'mumIsoPt_dr05/D'   )   
ntuple.Branch('mupIsoPt_dr05'    , mupIsoPt_dr05    , 'mupIsoPt_dr05/D'   )   
ntuple.Branch('kstTrkmIsoPt_dr05', kstTrkmIsoPt_dr05, 'kstTrkmIsoPt_dr05/D'   )   
ntuple.Branch('kstTrkpIsoPt_dr05', kstTrkpIsoPt_dr05, 'kstTrkpIsoPt_dr05/D'   )   

ntuple.Branch('charge_trig_matched', charge_trig_matched, 'charge_trig_matched/D')


if not skim:
  ntuple.Branch('mumNPixHits',                  mumNPixHits,                                  'mumNPixHits/D')
  ntuple.Branch('mumNTrkHits',                  mumNTrkHits,                                  'mumNTrkHits/D')
  ntuple.Branch('mupNPixHits',                  mupNPixHits,                                  'mupNPixHits/D')
  ntuple.Branch('mupNTrkHits',                  mupNTrkHits,                                  'mupNTrkHits/D')
  ntuple.Branch('mumKinkChi2',                  mumKinkChi2,                                  'mumKinkChi2/D')
  ntuple.Branch('mumFracHits',                  mumFracHits,                                  'mumFracHits/D')
  ntuple.Branch('mupKinkChi2',                  mupKinkChi2,                                  'mupKinkChi2/D')
  ntuple.Branch('mupFracHits',                  mupFracHits,                                  'mupFracHits/D')
  ntuple.Branch('mumMinIP',                     mumMinIP,                                     'mumMinIP/D')
  ntuple.Branch('mumMinIPS',                    mumMinIPS,                                    'mumMinIPS/D')
  ntuple.Branch('mupMinIP',                     mupMinIP,                                     'mupMinIP/D')
  ntuple.Branch('mupMinIPS',                    mupMinIPS,                                    'mupMinIPS/D')
  ntuple.Branch('mumGlobalMuonPromptTight',     mumGlobalMuonPromptTight      ,               'mumGlobalMuonPromptTight/D')
  ntuple.Branch('mumTMLastStationTight',        mumTMLastStationTight         ,               'mumTMLastStationTight/D')
  ntuple.Branch('mumTMLastStationLoose',        mumTMLastStationLoose         ,               'mumTMLastStationLoose/D')
  ntuple.Branch('mumTM2DCompatibilityTight',    mumTM2DCompatibilityTight     ,               'mumTM2DCompatibilityTight/D')
  ntuple.Branch('mumTM2DCompatibilityLoose',    mumTM2DCompatibilityLoose     ,               'mumTM2DCompatibilityLoose/D')
  ntuple.Branch('mumTMLastStationAngTight',     mumTMLastStationAngTight      ,               'mumTMLastStationAngTight/D')
  ntuple.Branch('mumTMLastStationAngLoose',     mumTMLastStationAngLoose      ,               'mumTMLastStationAngLoose/D')
  ntuple.Branch('mumTMOneStationAngTight',      mumTMOneStationAngTight       ,               'mumTMOneStationAngTight/D')
  ntuple.Branch('mumTMOneStationAngLoose',      mumTMOneStationAngLoose       ,               'mumTMOneStationAngLoose/D')
  ntuple.Branch('mumTrackerMuonArbitrated',     mumTrackerMuonArbitrated      ,               'mumTrackerMuonArbitrated/D')
  ntuple.Branch('mupGlobalMuonPromptTight',     mupGlobalMuonPromptTight      ,               'mupGlobalMuonPromptTight/D')
  ntuple.Branch('mupTMLastStationTight',        mupTMLastStationTight         ,               'mupTMLastStationTight/D')
  ntuple.Branch('mupTMLastStationLoose',        mupTMLastStationLoose         ,               'mupTMLastStationLoose/D')
  ntuple.Branch('mupTM2DCompatibilityTight',    mupTM2DCompatibilityTight     ,               'mupTM2DCompatibilityTight/D')
  ntuple.Branch('mupTM2DCompatibilityLoose',    mupTM2DCompatibilityLoose     ,               'mupTM2DCompatibilityLoose/D')
  ntuple.Branch('mupTMLastStationAngTight',     mupTMLastStationAngTight      ,               'mupTMLastStationAngTight/D')
  ntuple.Branch('mupTMLastStationAngLoose',     mupTMLastStationAngLoose      ,               'mupTMLastStationAngLoose/D')
  ntuple.Branch('mupTMOneStationAngTight',      mupTMOneStationAngTight       ,               'mupTMOneStationAngTight/D')
  ntuple.Branch('mupTMOneStationAngLoose',      mupTMOneStationAngLoose       ,               'mupTMOneStationAngLoose/D')
  ntuple.Branch('mupTrackerMuonArbitrated',     mupTrackerMuonArbitrated      ,               'mupTrackerMuonArbitrated/D')
  ntuple.Branch('kstTrkmTM2DCompatibilityTight',kstTrkmTM2DCompatibilityTight ,               'kstTrkmTM2DCompatibilityTight/D')
  ntuple.Branch('kstTrkmTM2DCompatibilityLoose',kstTrkmTM2DCompatibilityLoose ,               'kstTrkmTM2DCompatibilityLoose/D')
  ntuple.Branch('kstTrkpTM2DCompatibilityTight',kstTrkpTM2DCompatibilityTight ,               'kstTrkpTM2DCompatibilityTight/D')
  ntuple.Branch('kstTrkpTM2DCompatibilityLoose',kstTrkpTM2DCompatibilityLoose ,               'kstTrkpTM2DCompatibilityLoose/D')
  ntuple.Branch('kstTrkmTMLastStationAngTight', kstTrkmTMLastStationAngTight  ,               'kstTrkmTMLastStationAngTight/D')
  ntuple.Branch('kstTrkmTMLastStationAngLoose', kstTrkmTMLastStationAngLoose  ,               'kstTrkmTMLastStationAngLoose/D')
  ntuple.Branch('kstTrkpTMOneStationAngTight',  kstTrkpTMOneStationAngTight   ,               'kstTrkpTMOneStationAngTight/D')
  ntuple.Branch('kstTrkpTMOneStationAngLoose',  kstTrkpTMOneStationAngLoose   ,               'kstTrkpTMOneStationAngLoose/D')
  ntuple.Branch('kstTrkmTMOneStationAngTight',  kstTrkmTMOneStationAngTight   ,               'kstTrkmTMOneStationAngTight/D')
  ntuple.Branch('kstTrkmTMOneStationAngLoose',  kstTrkmTMOneStationAngLoose   ,               'kstTrkmTMOneStationAngLoose/D')
  ntuple.Branch('kstTrkpTMLastStationAngTight', kstTrkpTMLastStationAngTight  ,               'kstTrkpTMLastStationAngTight/D')
  ntuple.Branch('kstTrkpTMLastStationAngLoose', kstTrkpTMLastStationAngLoose  ,               'kstTrkpTMLastStationAngLoose/D')
  ntuple.Branch('kstTrkmTMLastStationTight',    kstTrkmTMLastStationTight     ,               'kstTrkmTMLastStationTight/D')
  ntuple.Branch('kstTrkmTMLastStationLoose',    kstTrkmTMLastStationLoose     ,               'kstTrkmTMLastStationLoose/D')
  ntuple.Branch('kstTrkpTMLastStationTight',    kstTrkpTMLastStationTight     ,               'kstTrkpTMLastStationTight/D')
  ntuple.Branch('kstTrkpTMLastStationLoose',    kstTrkpTMLastStationLoose     ,               'kstTrkpTMLastStationLoose/D')
  ntuple.Branch('kstTrkpMinIP',                 kstTrkpMinIP,                                 'kstTrkpMinIP/D')
  ntuple.Branch('kstTrkpMinIPS',                kstTrkpMinIPS,                                'kstTrkpMinIPS/D')
  ntuple.Branch('kstTrkmMinIP',                 kstTrkmMinIP,                                 'kstTrkmMinIP/D')
  ntuple.Branch('kstTrkmMinIPS',                kstTrkmMinIPS,                                'kstTrkmMinIPS/D')

  ntuple.Branch('mumIso'    ,  mumIso    )
  ntuple.Branch('mupIso'    ,  mupIso    )
  ntuple.Branch('kstTrkmIso',  kstTrkmIso)
  ntuple.Branch('kstTrkpIso',  kstTrkpIso)



numEvents = tree_lmnr.GetEntries()
print 'total number of events in tree:', numEvents

ROOT.gROOT.LoadMacro('FindValueFromVectorOfBool.h+')

progressbarWidth = 40
sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
sys.stdout.flush()                          # this forces to print the stdout buffer
sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['

tree_lmnr.SetBranchAddress('eventN',eventN)

for i, ev in enumerate(tree_lmnr):

    if i%int(numEvents/(progressbarWidth-1))==0:
        sys.stdout.write('+')
        sys.stdout.flush()

    for var in evbr:
        var[0] = 0
    
    if not len(ev.bMass) > 0:
        continue
        
    if not any( path in ev.TrigTable[0] for path in paths):
        continue     

    hlt_mums  = [ihlt for ihlt in ev.hltObjs if ihlt.pdgId == 13 ]
    hlt_mups  = [ihlt for ihlt in ev.hltObjs if ihlt.pdgId ==-13 ]
    hlt_trks  = [ihlt for ihlt in ev.hltObjs if abs(ihlt.pdgId) == 321]
    
    ## make all possible mumutk triplets from hlt
    triplets = list(itertools.product(hlt_mums,hlt_mups,hlt_trks))

    ## now loop on candidates per event
    for icand in range(len(ev.bMass)):
    
        for var in newbr:
            var[0] = -99.
        for var in isobr:
            var[0] = 0.
    
        if skimSoftMu and not (ev.mumNTrkLayers[icand] >= 6        and ev.mupNTrkLayers[icand] >= 6  and \
                               ev.mumNPixLayers[icand] >= 1        and ev.mupNPixLayers[icand] >= 1  and \
                               ev.mumdxyBS[icand] < 0.3            and ev.mupdxyBS[icand] < 0.3      and \
                               ev.mumdzBS[icand] < 20              and ev.mupdzBS[icand]  < 20       and \
                               ROOT.FindValueFromVectorOfBool(ev.mumHighPurity, icand) == 1          and \
                               ROOT.FindValueFromVectorOfBool(ev.mupHighPurity, icand) == 1        ):
            continue


        ## trigger match: both muons should be matched + one track
        ## save the charge of the matched track to then apply pT cuts only to one of them
        charge_matched = findTriggerMatching(triplets,
                                             ev.rawmumEta[icand],     ev.rawmumPhi[icand],     ev.rawmumPt[icand],
                                             ev.rawmupEta[icand],     ev.rawmupPhi[icand],     ev.rawmupPt[icand],
                                             ev.rawkstTrkmEta[icand], ev.rawkstTrkmPhi[icand], ev.rawkstTrkmPt[icand],
                                             ev.rawkstTrkpEta[icand], ev.rawkstTrkpPhi[icand], ev.rawkstTrkpPt[icand]
                                             ) 
        if charge_matched == 0:
            continue
   

    

        ## per event quantities
        runN[0]                        = ev.runN
        recoVtxN[0]                    = ev.recoVtxN
        trig[0]                        = paths.index(ev.TrigTable[0].split('_v')[0])

        l1_00_1p5[0]                   = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4')
        l1_00_1p4[0]                   = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4')
        l1_4_4[0]                      = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4')
        l1_4p5_4p5[0]                  = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2')

#         l1_11_4[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_11_4')
#         l1_12_5[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_12_5')
#         l1_10_0[0]                     = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu_10_0_dEta_Max1p8')
#         l1_00[0]                       = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p6_dEta_Max1p8')
#         l1_00_OS[0]                    = findFiringL1(ev.L1Table, ev.L1Prescales, 'L1_DoubleMu0er1p6_dEta_Max1p8_OS')

        bsX[0]                         = ev.bsX
        bsY[0]                         = ev.bsY

        bPt[0]                         = computePt( ev.bPx      [icand], ev.bPy       [icand] ) 
        kstPt[0]                       = computePt( ev.kstPx    [icand], ev.kstPy     [icand] ) 
        mumuPt[0]                      = computePt( ev.mumuPx   [icand], ev.mumuPy    [icand] ) 
        mumPt[0]                       = computePt( ev.mumPx    [icand], ev.mumPy     [icand] ) 
        mupPt[0]                       = computePt( ev.mupPx    [icand], ev.mupPy     [icand] ) 
        kstTrkmPt[0]                   = computePt( ev.kstTrkmPx[icand], ev.kstTrkmPy [icand] ) 
        kstTrkpPt[0]                   = computePt( ev.kstTrkpPx[icand], ev.kstTrkpPy [icand] ) 

        bPhi[0]                        = computePhi( ev.bPx      [icand], ev.bPy       [icand] )
        kstPhi[0]                      = computePhi( ev.kstPx    [icand], ev.kstPy     [icand] )
        mumuPhi[0]                     = computePhi( ev.mumuPx   [icand], ev.mumuPy    [icand] )
        mumPhi[0]                      = computePhi( ev.mumPx    [icand], ev.mumPy     [icand] )
        mupPhi[0]                      = computePhi( ev.mupPx    [icand], ev.mupPy     [icand] )
        kstTrkmPhi[0]                  = computePhi( ev.kstTrkmPx[icand], ev.kstTrkmPy [icand] )
        kstTrkpPhi[0]                  = computePhi( ev.kstTrkpPx[icand], ev.kstTrkpPy [icand] )

        bEta[0]                        = computeEta( ev.bPx      [icand], ev.bPy       [icand], ev.bPz       [icand] )
        kstEta[0]                      = computeEta( ev.kstPx    [icand], ev.kstPy     [icand], ev.kstPz     [icand] )
        mumuEta[0]                     = computeEta( ev.mumuPx   [icand], ev.mumuPy    [icand], ev.mumuPz    [icand] )
        mumEta[0]                      = computeEta( ev.mumPx    [icand], ev.mumPy     [icand], ev.mumPz     [icand] )
        mupEta[0]                      = computeEta( ev.mupPx    [icand], ev.mupPy     [icand], ev.mupPz     [icand] )
        kstTrkmEta[0]                  = computeEta( ev.kstTrkmPx[icand], ev.kstTrkmPy [icand], ev.kstTrkmPz [icand] )
        kstTrkpEta[0]                  = computeEta( ev.kstTrkpPx[icand], ev.kstTrkpPy [icand], ev.kstTrkpPz [icand] )

        ## per-candidate quantities
        bMass[0]                       = ev.bMass[icand]
        bMassE[0]                      = ev.bMassE[icand]
        bBarMass[0]                    = ev.bBarMass[icand]
        bBarMassE[0]                   = ev.bBarMassE[icand]
        bVtxCL[0]                      = ev.bVtxCL[icand]
        bVtxX[0]                       = ev.bVtxX[icand]
        bVtxY[0]                       = ev.bVtxY[icand]
        bVtxZ[0]                       = ev.bVtxZ[icand]
        
        bCosAlphaBS[0]                 = ev.bCosAlphaBS[icand]
        bCosAlphaBSE[0]                = ev.bCosAlphaBSE[icand]
        bLBS[0]                        = ev.bLBS[icand]
        bLBSE[0]                       = ev.bLBSE[icand]
        bDCABS[0]                      = ev.bDCABS[icand]
        bDCABSE[0]                     = ev.bDCABSE[icand]
        kstMass[0]                     = ev.kstMass[icand]
        kstMassE[0]                    = ev.kstMassE[icand]
        kstBarMass[0]                  = ev.kstBarMass[icand]
        kstBarMassE[0]                 = ev.kstBarMassE[icand]
        kstVtxCL[0]                    = ev.kstVtxCL[icand]
        kstVtxX[0]                     = ev.kstVtxX[icand]
        kstVtxY[0]                     = ev.kstVtxY[icand]
        kstVtxZ[0]                     = ev.kstVtxZ[icand]
        kkMass[0]                      = computeInvMass(ev.kstTrkmPx[icand], ev.kstTrkmPy[icand], ev.kstTrkmPz[icand], kaonMass,
                                                        ev.kstTrkpPx[icand], ev.kstTrkpPy[icand], ev.kstTrkpPz[icand], kaonMass )

        
        mumuMass[0]                    = ev.mumuMass[icand]
        mumuMassE[0]                   = ev.mumuMassE[icand]
        mumuVtxCL[0]                   = ev.mumuVtxCL[icand]
        mumuVtxX[0]                    = ev.mumuVtxX[icand]
        mumuVtxY[0]                    = ev.mumuVtxY[icand]
        mumuVtxZ[0]                    = ev.mumuVtxZ[icand]
        mumuCosAlphaBS[0]              = ev.mumuCosAlphaBS[icand]
        mumuCosAlphaBSE[0]             = ev.mumuCosAlphaBSE[icand]
        mumuLBS[0]                     = ev.mumuLBS[icand]
        mumuLBSE[0]                    = ev.mumuLBSE[icand]
        mumuDCA[0]                     = ev.mumuDCA[icand]
        
        mumHighPurity[0]               = ROOT.FindValueFromVectorOfBool(ev.mumHighPurity, icand) 
        mumCL[0]                       = ev.mumCL[icand]
        mumNormChi2[0]                 = ev.mumNormChi2[icand]
        mumDCABS[0]                    = ev.mumDCABS[icand]
        mumDCABSE[0]                   = ev.mumDCABSE[icand]
        mumdxyBS[0]                    = ev.mumdxyBS[icand]
        mumdzBS[0]                     = ev.mumdzBS[icand]
        mumMinIP2D[0]                  = ev.mumMinIP2D[icand]
        mumMinIP2DE[0]                 = ev.mumMinIP2DE[icand]
        mumNPixLayers[0]               = ev.mumNPixLayers[icand]
        mumNTrkLayers[0]               = ev.mumNTrkLayers[icand]
        mumNMuonHits[0]                = ev.mumNMuonHits[icand]
        mumNMatchStation[0]            = ev.mumNMatchStation[icand]
        
        mumCategoryDict = muonCategory(ev.mumCat[icand])
        mumGlobalMuon            [0] = mumCategoryDict['GlobalMuon']
        mumTrackerMuon           [0] = mumCategoryDict['TrackerMuon']
        mumStandAloneMuon        [0] = mumCategoryDict['StandAloneMuon']
        mumTMOneStationTight     [0] = mumCategoryDict['TMOneStationTight']
        mumTMOneStationLoose     [0] = mumCategoryDict['TMOneStationLoose']        

        mupHighPurity[0]               = ROOT.FindValueFromVectorOfBool(ev.mupHighPurity, icand)
        mupCL[0]                       = ev.mupCL[icand]
        mupNormChi2[0]                 = ev.mupNormChi2[icand]
        mupDCABS[0]                    = ev.mupDCABS[icand]
        mupDCABSE[0]                   = ev.mupDCABSE[icand]
        mupdxyBS[0]                    = ev.mupdxyBS[icand]
        mupdzBS[0]                     = ev.mupdzBS[icand]
        mupMinIP2D[0]                  = ev.mupMinIP2D[icand]
        mupMinIP2DE[0]                 = ev.mupMinIP2DE[icand]
        mupNPixLayers[0]               = ev.mupNPixLayers[icand]
        mupNTrkLayers[0]               = ev.mupNTrkLayers[icand]
        mupNMuonHits[0]                = ev.mupNMuonHits[icand]
        mupNMatchStation[0]            = ev.mupNMatchStation[icand]
        
        mupCategoryDict = muonCategory(ev.mupCat[icand])
        mupGlobalMuon            [0] = mupCategoryDict['GlobalMuon']
        mupTrackerMuon           [0] = mupCategoryDict['TrackerMuon']
        mupStandAloneMuon        [0] = mupCategoryDict['StandAloneMuon']
        mupTMOneStationTight     [0] = mupCategoryDict['TMOneStationTight']
        mupTMOneStationLoose     [0] = mupCategoryDict['TMOneStationLoose']
        
        
        kstTrkmHighPurity[0]           = ROOT.FindValueFromVectorOfBool(ev.kstTrkmHighPurity, icand)
        kstTrkmCL[0]                   = ev.kstTrkmCL[icand]
        kstTrkmNormChi2[0]             = ev.kstTrkmNormChi2[icand]
        kstTrkmPz[0]                   = ev.kstTrkmPz[icand]
        kstTrkmDCABS[0]                = ev.kstTrkmDCABS[icand]     ## theDCAXBS.perigeeParameters().transverseImpactParameter();
        kstTrkmDCABSE[0]               = ev.kstTrkmDCABSE[icand]
        kstTrkmFracHits[0]             = ev.kstTrkmFracHits[icand]
        kstTrkmNPixHits[0]             = ev.kstTrkmNPixHits[icand]
        kstTrkmNPixLayers[0]           = ev.kstTrkmNPixLayers[icand]
        kstTrkmNTrkHits[0]             = ev.kstTrkmNTrkHits[icand]
        kstTrkmNTrkLayers[0]           = ev.kstTrkmNTrkLayers[icand]
        kstTrkmMinIP2D[0]              = ev.kstTrkmMinIP2D[icand]
        kstTrkmMinIP2DE[0]             = ev.kstTrkmMinIP2DE[icand]

        trkmCategoryDict = muonCategory(ev.kstTrkmMuMatch[icand])
        kstTrkmGlobalMuon            [0] = trkmCategoryDict['GlobalMuon']
        kstTrkmTrackerMuon           [0] = trkmCategoryDict['TrackerMuon']
        kstTrkmStandAloneMuon        [0] = trkmCategoryDict['StandAloneMuon']
        kstTrkmTrackerMuonArbitrated [0] = trkmCategoryDict['TrackerMuonArbitrated']
        kstTrkmTMOneStationTight     [0] = trkmCategoryDict['TMOneStationTight']
        kstTrkmTMOneStationLoose     [0] = trkmCategoryDict['TMOneStationLoose']

        kstTrkpHighPurity[0]           = ROOT.FindValueFromVectorOfBool(ev.kstTrkpHighPurity, icand)
        kstTrkpCL[0]                   = ev.kstTrkpCL[icand]
        kstTrkpNormChi2[0]             = ev.kstTrkpNormChi2[icand]
        kstTrkpPz[0]                   = ev.kstTrkpPz[icand]
        kstTrkpDCABS[0]                = ev.kstTrkpDCABS[icand]
        kstTrkpDCABSE[0]               = ev.kstTrkpDCABSE[icand]
        kstTrkpFracHits[0]             = ev.kstTrkpFracHits[icand]
        kstTrkpNPixHits[0]             = ev.kstTrkpNPixHits[icand]
        kstTrkpNPixLayers[0]           = ev.kstTrkpNPixLayers[icand]
        kstTrkpNTrkHits[0]             = ev.kstTrkpNTrkHits[icand]
        kstTrkpNTrkLayers[0]           = ev.kstTrkpNTrkLayers[icand]
        kstTrkpMinIP2D[0]              = ev.kstTrkpMinIP2D[icand]
        kstTrkpMinIP2DE[0]             = ev.kstTrkpMinIP2DE[icand]

        trkpCategoryDict = muonCategory(ev.kstTrkpMuMatch[icand])
        kstTrkpGlobalMuon            [0] = trkpCategoryDict['GlobalMuon']
        kstTrkpTrackerMuon           [0] = trkpCategoryDict['TrackerMuon']
        kstTrkpStandAloneMuon        [0] = trkpCategoryDict['StandAloneMuon']
        kstTrkpTrackerMuonArbitrated [0] = trkpCategoryDict['TrackerMuonArbitrated']
        kstTrkpTMOneStationTight     [0] = trkpCategoryDict['TMOneStationTight']
        kstTrkpTMOneStationLoose     [0] = trkpCategoryDict['TMOneStationLoose']

        tagB0[0]                      = FlavorTagger(ev.kstMass[icand], ev.kstBarMass[icand])  ## 1 if B0, 0 if B0bar

        bMinusVtxCL[0]               = ev.bMinusVtxCL[icand]      
        bMinusCosAlphaBS[0]          = ev.bMinusCosAlphaBS[icand]
        bPlusVtxCL[0]                = ev.bPlusVtxCL[icand]    
        bPlusCosAlphaBS[0]           = ev.bPlusCosAlphaBS[icand]  

        charge_trig_matched[0]  =  charge_matched    

        cos_theta_l[0], cos_theta_k[0], phi_kst_mumu[0] = addVars(
        	tagB0[0],
        	ev.mumuPx[icand],     ev.mumuPy[icand],    ev.mumuPz[icand], ev.mumuMass[icand],
        	ev.mumPx[icand],      ev.mumPy[icand],     ev.mumPz[icand],  
        	ev.mupPx[icand],      ev.mupPy[icand],     ev.mupPz[icand],  
        	ev.kstTrkmPx[icand],  ev.kstTrkmPy[icand], ev.kstTrkmPz[icand],
        	ev.kstTrkpPx[icand],  ev.kstTrkpPy[icand], ev.kstTrkpPz[icand],
        	ev.kstPx[icand],      ev.kstPy[icand],     ev.kstPz[icand],   ev.kstMass[icand], ev.kstBarMass[icand],
        	ev.bPx[icand],        ev.bPy[icand],       ev.bPz[icand],     ev.bMass[icand],   ev.bBarMass[icand]
        )

        mmk1[0], mmk2[0] = addMMKVars(
                          mumPt[0],      mumEta[0],      mumPhi[0],  
                          mupPt[0],      mupEta[0],      mupPhi[0],  
                          kstTrkpPt[0],  kstTrkpEta[0],  kstTrkpPhi[0],
                          kstTrkmPt[0],  kstTrkmEta[0],  kstTrkmPhi[0]
                          )
        dR_mum_trkm[0] = addDR(mumEta[0], mumPhi[0], kstTrkmEta[0], kstTrkmPhi[0])
        dR_mup_trkp[0] = addDR(mupEta[0], mupPhi[0], kstTrkpEta[0], kstTrkpPhi[0])



        ###########  mu - isolation ####################
        val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
        val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
        val_isoPt_dr05 = 0;     val_isoP_dr05  = 0;
        
        if len( ev.mumIsoPt[icand] ) > 0:
            for j,isoi in enumerate(ev.mumIsoPt[icand]):
                if ev.mumIsodR[icand][j] < 0.3:
                    val_isoPt_dr03 += isoi
                if ev.mumIsodR[icand][j] < 0.4:
                    val_isoPt_dr04 += isoi
                if ev.mumIsodR[icand][j] < 0.5:
                    val_isoPt_dr05 += isoi

            mumIsoPt_dr03[0] = val_isoPt_dr03
            mumIsoPt_dr04[0] = val_isoPt_dr04
            mumIsoPt_dr05[0] = val_isoPt_dr05

        ###########  mu + isolation ####################
        val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
        val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
        val_isoPt_dr05 = 0;     val_isoP_dr05  = 0;
        
        if len( ev.mupIsoPt[icand] ) > 0:
            for j,isoi in enumerate(ev.mupIsoPt[icand]):
                if ev.mupIsodR[icand][j] < 0.3:
                    val_isoPt_dr03 += isoi
                if ev.mupIsodR[icand][j] < 0.4:
                    val_isoPt_dr04 += isoi
                if ev.mupIsodR[icand][j] < 0.5:
                    val_isoPt_dr05 += isoi

            mupIsoPt_dr03[0] = val_isoPt_dr03
            mupIsoPt_dr04[0] = val_isoPt_dr04
            mupIsoPt_dr05[0] = val_isoPt_dr05

 
        ###########  trk - isolation ####################
        val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
        val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
        val_isoPt_dr05 = 0;     val_isoP_dr05  = 0;
        
        if len( ev.kstTrkmIsoPt[icand] ) > 0:
            for j,isoi in enumerate(ev.kstTrkmIsoPt[icand]):
                if ev.kstTrkmIsodR[icand][j] < 0.3:
                    val_isoPt_dr03 += isoi
                if ev.kstTrkmIsodR[icand][j] < 0.4:
                    val_isoPt_dr04 += isoi
                if ev.kstTrkmIsodR[icand][j] < 0.5:
                    val_isoPt_dr05 += isoi

            kstTrkmIsoPt_dr03[0] = val_isoPt_dr03
            kstTrkmIsoPt_dr04[0] = val_isoPt_dr04
            kstTrkmIsoPt_dr05[0] = val_isoPt_dr05

        ###########  trk + isolation ####################
        val_isoPt_dr03 = 0;     val_isoP_dr03  = 0;
        val_isoPt_dr04 = 0;     val_isoP_dr04  = 0;
        val_isoPt_dr05 = 0;     val_isoP_dr05  = 0;

        if len( ev.kstTrkpIsoPt[icand] ) > 0:
            for j,isoi in enumerate(ev.kstTrkpIsoPt[icand]):
                if ev.kstTrkpIsodR[icand][j] < 0.3:
                    val_isoPt_dr03 += isoi
                if ev.kstTrkpIsodR[icand][j] < 0.4:
                    val_isoPt_dr04 += isoi
                if ev.kstTrkpIsodR[icand][j] < 0.5:
                    val_isoPt_dr05 += isoi

            kstTrkpIsoPt_dr03[0] = val_isoPt_dr03
            kstTrkpIsoPt_dr04[0] = val_isoPt_dr04
            kstTrkpIsoPt_dr05[0] = val_isoPt_dr05


        if not skim:
            mumNPixHits[0]                 = ev.mumNPixHits[icand]
            mumNTrkHits[0]                 = ev.mumNTrkHits[icand]
            mumKinkChi2[0]                 = ev.mumKinkChi2[icand]
            mumFracHits[0]                 = ev.mumFracHits[icand]
            mumMinIP[0]                    = ev.mumMinIP[icand]
            mumMinIPS[0]                   = ev.mumMinIPS[icand]
            mumTMLastStationTight    [0]   = mumCategoryDict['TMLastStationTight']
            mumTMLastStationLoose    [0]   = mumCategoryDict['TMLastStationLoose']
            mumTrackerMuonArbitrated [0]   = mumCategoryDict['TrackerMuonArbitrated']
            mumGlobalMuonPromptTight [0]   = mumCategoryDict['GlobalMuonPromptTight']
            mumTM2DCompatibilityTight[0]   = mumCategoryDict['TM2DCompatibilityTight']
            mumTM2DCompatibilityLoose[0]   = mumCategoryDict['TM2DCompatibilityLoose']
            mumTMLastStationAngTight [0]   = mumCategoryDict['TMLastStationAngTight']
            mumTMLastStationAngLoose [0]   = mumCategoryDict['TMLastStationAngLoose']
            mumTMOneStationAngTight  [0]   = mumCategoryDict['TMOneStationAngTight']
            mumTMOneStationAngLoose  [0]   = mumCategoryDict['TMOneStationAngLoose']
            mupKinkChi2[0]                 = ev.mupKinkChi2[icand]
            mupFracHits[0]                 = ev.mupFracHits[icand]
            mupGlobalMuonPromptTight [0]   = mupCategoryDict['GlobalMuonPromptTight']
            mupTMLastStationTight    [0]   = mupCategoryDict['TMLastStationTight']
            mupTMLastStationLoose    [0]   = mupCategoryDict['TMLastStationLoose']
            mupTM2DCompatibilityTight[0]   = mupCategoryDict['TM2DCompatibilityTight']
            mupTM2DCompatibilityLoose[0]   = mupCategoryDict['TM2DCompatibilityLoose']
            mupTMLastStationAngTight [0]   = mupCategoryDict['TMLastStationAngTight']
            mupTMLastStationAngLoose [0]   = mupCategoryDict['TMLastStationAngLoose']
            mupTMOneStationAngTight  [0]   = mupCategoryDict['TMOneStationAngTight']
            mupTMOneStationAngLoose  [0]   = mupCategoryDict['TMOneStationAngLoose']
            mupTrackerMuonArbitrated [0] = mupCategoryDict['TrackerMuonArbitrated']
            mupNPixHits[0]                 = ev.mupNPixHits[icand]
            mupNTrkHits[0]                 = ev.mupNTrkHits[icand]
            mupMinIP[0]                    = ev.mupMinIP[icand]
            mupMinIPS[0]                   = ev.mupMinIPS[icand]

            kstTrkmGlobalMuonPromptTight [0] = trkmCategoryDict['GlobalMuonPromptTight']
            kstTrkmTMLastStationTight    [0] = trkmCategoryDict['TMLastStationTight']
            kstTrkmTMLastStationLoose    [0] = trkmCategoryDict['TMLastStationLoose']
            kstTrkmTM2DCompatibilityTight[0] = trkmCategoryDict['TM2DCompatibilityTight']
            kstTrkmTM2DCompatibilityLoose[0] = trkmCategoryDict['TM2DCompatibilityLoose']
            kstTrkmTMLastStationAngTight [0] = trkmCategoryDict['TMLastStationAngTight']
            kstTrkmTMLastStationAngLoose [0] = trkmCategoryDict['TMLastStationAngLoose']
            kstTrkmTMOneStationAngTight  [0] = trkmCategoryDict['TMOneStationAngTight']
            kstTrkmTMOneStationAngLoose  [0] = trkmCategoryDict['TMOneStationAngLoose']
  
            kstTrkmMinIP[0]                = ev.kstTrkmMinIP[icand]
            kstTrkmMinIPS[0]               = ev.kstTrkmMinIPS[icand]
            kstTrkpMinIP[0]                = ev.kstTrkpMinIP[icand]
            kstTrkpMinIPS[0]               = ev.kstTrkpMinIPS[icand]
            kstTrkpGlobalMuonPromptTight [0] = trkpCategoryDict['GlobalMuonPromptTight']
            kstTrkpTMLastStationTight    [0] = trkpCategoryDict['TMLastStationTight']
            kstTrkpTMLastStationLoose    [0] = trkpCategoryDict['TMLastStationLoose']
            kstTrkpTM2DCompatibilityTight[0] = trkpCategoryDict['TM2DCompatibilityTight']
            kstTrkpTM2DCompatibilityLoose[0] = trkpCategoryDict['TM2DCompatibilityLoose']
            kstTrkpTMLastStationAngTight [0] = trkpCategoryDict['TMLastStationAngTight']
            kstTrkpTMLastStationAngLoose [0] = trkpCategoryDict['TMLastStationAngLoose']
            kstTrkpTMOneStationAngTight  [0] = trkpCategoryDict['TMOneStationAngTight']
            kstTrkpTMOneStationAngLoose  [0] = trkpCategoryDict['TMOneStationAngLoose']

            mumIso    .clear()
            mupIso    .clear()
            kstTrkmIso.clear()
            kstTrkpIso.clear()
    
            if len( ev.mumIso[icand] ) > 0:
                for i,isoi in enumerate(ev.mumIso[icand]):
                    mumIso.push_back(isoi)
            else:
                mumIso.push_back(0)        
    
            if len( ev.mupIso[icand] ) > 0:
                for i,isoi in enumerate(ev.mupIso[icand]):
                    mupIso.push_back(isoi)
            else:
                mupIso.push_back(0)        
    
            if len( ev.kstTrkmIso[icand] ) > 0:
                for i,isoi in enumerate(ev.kstTrkmIso[icand]):
                    kstTrkmIso.push_back(isoi)
            else:
                kstTrkmIso.push_back(0)        
    
            if len( ev.kstTrkpIso[icand] ) > 0:
                for i,isoi in enumerate(ev.kstTrkpIso[icand]):
                    kstTrkpIso.push_back(isoi)
            else:
                kstTrkpIso.push_back(0)        
    

        ntuple. Fill()


sys.stdout.write('\n')

file_out.cd()
ntuple.Write()
file_out.Close()
    
