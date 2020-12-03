from argparse import ArgumentParser
from samples import *

parser = ArgumentParser()
parser.add_argument("sample",                                       help = "sample", nargs = '+', choices = samples.keys(), default = 'BJpsiK_ee_mc_2019Oct25' )
parser.add_argument("-n"  , "--number"     , dest = "number"     ,  help = "number of file"        , default = '1')
parser.add_argument("-f"  , "--folder"     , dest = "folder"     ,  help = "name of output folder" , default = 'default_folder')
parser.add_argument("-g"  , "--dogen"      , dest = "dogen"      ,  help = "produce gen ntuples "  , default=False, action='store_true'      )
parser.add_argument("-c"  , "--chan"       , dest = "thechan"    ,  help = "LMNR, PSI"                   , default = 'LMNR' ) ## not used here


args = parser.parse_args()
if not args.number:   
  parser.error('Number filename not given')

import ROOT
from DataFormats.FWLite import Events
from array import array
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
import numpy as np

import sys, math, itertools
sys.path.append('/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu')
from utils.angular_vars import *
from utils.utils import *


type       = args.thechan  ## LMNR, JPsi
tree_lmnr = ROOT.TChain('B0KstMuMu/B0KstMuMuNTuple')

indir = '0000'
if int(args.number) >= 1000:  indir = '0001'
if int(args.number) >= 2000:  indir = '0002'
if int(args.number) >= 3000:  indir = '0003'

tree_lmnr.Add('%s/%s/B0ToKstMuMu_%s.root'        %( samples[args.sample[0]]['path'],indir,args.number))


file_out = ROOT.TFile('gen_ntuple_%s_noFilter_%s.root'%(type,args.number), 'recreate')
ntuple   = ROOT.TTree( 'ntuple', 'ntuple' )


print '@@@@@@@@@@@     flattening ntuple     @@@@@@@@@@@'
print 'input file: ', tree_lmnr.GetFile()

try:
  tree_lmnr.GetEntries()
except:
  print 'ntuple tree not found'
  exit()


genbr  =  []

'''
for i in branches:
    print "{NAME}\t = array('d', [-99.]);  recobr.append({NAME})".format(NAME=i).expandtabs(40)
'''
eventN                           = array('L', [99]);  
runN                             = array('d', [-99.]);  genbr.append(runN)
genSignal                        = array('d', [-99.]);  genbr.append(genSignal)
genSignHasFSR                    = array('d', [-99.]);  genbr.append(genSignHasFSR)
genSignKstHasFSR                 = array('d', [-99.]);  genbr.append(genSignKstHasFSR)
genSignPsiHasFSR                 = array('d', [-99.]);  genbr.append(genSignPsiHasFSR)
# genPriVtxX                       = array('d', [-99.]);  genbr.append(genPriVtxX)
# genPriVtxY                       = array('d', [-99.]);  genbr.append(genPriVtxY)
# genPriVtxZ                       = array('d', [-99.]);  genbr.append(genPriVtxZ)
genB0Mass                        = array('d', [-99.]);  genbr.append(genB0Mass)
# genB0Px                          = array('d', [-99.]);  genbr.append(genB0Px)
# genB0Py                          = array('d', [-99.]);  genbr.append(genB0Py)
# genB0Pz                          = array('d', [-99.]);  genbr.append(genB0Pz)
# genB0VtxX                        = array('d', [-99.]);  genbr.append(genB0VtxX)
# genB0VtxY                        = array('d', [-99.]);  genbr.append(genB0VtxY)
# genB0VtxZ                        = array('d', [-99.]);  genbr.append(genB0VtxZ)
genKstMass                       = array('d', [-99.]);  genbr.append(genKstMass)
# genKstPx                         = array('d', [-99.]);  genbr.append(genKstPx)
# genKstPy                         = array('d', [-99.]);  genbr.append(genKstPy)
# genKstPz                         = array('d', [-99.]);  genbr.append(genKstPz)
# genKstVtxX                       = array('d', [-99.]);  genbr.append(genKstVtxX)
# genKstVtxY                       = array('d', [-99.]);  genbr.append(genKstVtxY)
# genKstVtxZ                       = array('d', [-99.]);  genbr.append(genKstVtxZ)
genPsiMass                       = array('d', [-99.]);  genbr.append(genPsiMass)
# genPsiVtxX                       = array('d', [-99.]);  genbr.append(genPsiVtxX)
# genPsiVtxY                       = array('d', [-99.]);  genbr.append(genPsiVtxY)
# genPsiVtxZ                       = array('d', [-99.]);  genbr.append(genPsiVtxZ)
genMumMother                     = array('d', [-99.]);  genbr.append(genMumMother)
# genMumPx                         = array('d', [-99.]);  genbr.append(genMumPx)
# genMumPy                         = array('d', [-99.]);  genbr.append(genMumPy)
# genMumPz                         = array('d', [-99.]);  genbr.append(genMumPz)
genMupMother                     = array('d', [-99.]);  genbr.append(genMupMother)
# genMupPx                         = array('d', [-99.]);  genbr.append(genMupPx)
# genMupPy                         = array('d', [-99.]);  genbr.append(genMupPy)
# genMupPz                         = array('d', [-99.]);  genbr.append(genMupPz)
genKstTrkmMother                 = array('d', [-99.]);  genbr.append(genKstTrkmMother)
genKstTrkmID                     = array('d', [-99.]);  genbr.append(genKstTrkmID)
# genKstTrkmPx                     = array('d', [-99.]);  genbr.append(genKstTrkmPx)
# genKstTrkmPy                     = array('d', [-99.]);  genbr.append(genKstTrkmPy)
# genKstTrkmPz                     = array('d', [-99.]);  genbr.append(genKstTrkmPz)
genKstTrkpMother                 = array('d', [-99.]);  genbr.append(genKstTrkpMother)
genKstTrkpID                     = array('d', [-99.]);  genbr.append(genKstTrkpID)
# genKstTrkpPx                     = array('d', [-99.]);  genbr.append(genKstTrkpPx)
# genKstTrkpPy                     = array('d', [-99.]);  genbr.append(genKstTrkpPy)
# genKstTrkpPz                     = array('d', [-99.]);  genbr.append(genKstTrkpPz)

genbPt                              = array('d',[-99.]); genbr.append( genbPt          )   
genkstPt                            = array('d',[-99.]); genbr.append( genkstPt        )   
genmumuPt                           = array('d',[-99.]); genbr.append( genmumuPt       )   
genmumPt                            = array('d',[-99.]); genbr.append( genmumPt        )   
genmupPt                            = array('d',[-99.]); genbr.append( genmupPt        )   
genkstTrkmPt                        = array('d',[-99.]); genbr.append( genkstTrkmPt    )   
genkstTrkpPt                        = array('d',[-99.]); genbr.append( genkstTrkpPt    )   
genbPhi                             = array('d',[-99.]); genbr.append( genbPhi         )   
genkstPhi                           = array('d',[-99.]); genbr.append( genkstPhi       )   
genmumuPhi                          = array('d',[-99.]); genbr.append( genmumuPhi      )   
genmumPhi                           = array('d',[-99.]); genbr.append( genmumPhi       )   
genmupPhi                           = array('d',[-99.]); genbr.append( genmupPhi       )   
genkstTrkmPhi                       = array('d',[-99.]); genbr.append( genkstTrkmPhi   )   
genkstTrkpPhi                       = array('d',[-99.]); genbr.append( genkstTrkpPhi   )   
genbEta                             = array('d',[-99.]); genbr.append( genbEta         )   
genkstEta                           = array('d',[-99.]); genbr.append( genkstEta       )   
genmumuEta                          = array('d',[-99.]); genbr.append( genmumuEta      )   
genmumEta                           = array('d',[-99.]); genbr.append( genmumEta       )   
genmupEta                           = array('d',[-99.]); genbr.append( genmupEta       )   
genkstTrkmEta                       = array('d',[-99.]); genbr.append( genkstTrkmEta   )   
genkstTrkpEta                       = array('d',[-99.]); genbr.append( genkstTrkpEta   )   
genQ2                               = array('d', [-99.]);  genbr.append( genQ2           )   
genQ                                = array('d', [-99.]);  genbr.append( genQ            )   

cos_theta_l                     = array('d',[-99.]); genbr.append( cos_theta_l )   
cos_theta_k                     = array('d',[-99.]); genbr.append( cos_theta_k )   
phi_kst_mumu                    = array('d',[-99.]); genbr.append( phi_kst_mumu)   


'''
for i in branches:
    print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/D')".format(NAME=i).expandtabs(40)
'''
ntuple.Branch('runN',                   runN,                                   'runN/D')
ntuple.Branch('eventN',                 eventN,                                 'eventN/L')
ntuple.Branch('genSignal',              genSignal,                              'genSignal/D')
ntuple.Branch('genSignHasFSR',          genSignHasFSR,                          'genSignHasFSR/D')
ntuple.Branch('genSignKstHasFSR',       genSignKstHasFSR,                       'genSignKstHasFSR/D')
ntuple.Branch('genSignPsiHasFSR',       genSignPsiHasFSR,                       'genSignPsiHasFSR/D')
# ntuple.Branch('genPriVtxX',             genPriVtxX,                             'genPriVtxX/D')
# ntuple.Branch('genPriVtxY',             genPriVtxY,                             'genPriVtxY/D')
# ntuple.Branch('genPriVtxZ',             genPriVtxZ,                             'genPriVtxZ/D')
ntuple.Branch('genB0Mass',              genB0Mass,                              'genB0Mass/D')
# ntuple.Branch('genB0Px',                genB0Px,                                'genB0Px/D')
# ntuple.Branch('genB0Py',                genB0Py,                                'genB0Py/D')
# ntuple.Branch('genB0Pz',                genB0Pz,                                'genB0Pz/D')
# ntuple.Branch('genB0VtxX',              genB0VtxX,                              'genB0VtxX/D')
# ntuple.Branch('genB0VtxY',              genB0VtxY,                              'genB0VtxY/D')
# ntuple.Branch('genB0VtxZ',              genB0VtxZ,                              'genB0VtxZ/D')
ntuple.Branch('genKstMass',             genKstMass,                             'genKstMass/D')
# ntuple.Branch('genKstPx',               genKstPx,                               'genKstPx/D')
# ntuple.Branch('genKstPy',               genKstPy,                               'genKstPy/D')
# ntuple.Branch('genKstPz',               genKstPz,                               'genKstPz/D')
# ntuple.Branch('genKstVtxY',             genKstVtxY,                             'genKstVtxY/D')
ntuple.Branch('genPsiMass',             genPsiMass,                             'genPsiMass/D')
# ntuple.Branch('genPsiVtxX',             genPsiVtxX,                             'genPsiVtxX/D')
# ntuple.Branch('genPsiVtxY',             genPsiVtxY,                             'genPsiVtxY/D')
# ntuple.Branch('genPsiVtxZ',             genPsiVtxZ,                             'genPsiVtxZ/D')
ntuple.Branch('genMumMother',           genMumMother,                           'genMumMother/D')
# ntuple.Branch('genMumPx',               genMumPx,                               'genMumPx/D')
# ntuple.Branch('genMumPy',               genMumPy,                               'genMumPy/D')
# ntuple.Branch('genMumPz',               genMumPz,                               'genMumPz/D')
ntuple.Branch('genMupMother',           genMupMother,                           'genMupMother/D')
# ntuple.Branch('genMupPx',               genMupPx,                               'genMupPx/D')
# ntuple.Branch('genMupPy',               genMupPy,                               'genMupPy/D')
# ntuple.Branch('genMupPz',               genMupPz,                               'genMupPz/D')
ntuple.Branch('genKstTrkmMother',       genKstTrkmMother,                       'genKstTrkmMother/D')
# ntuple.Branch('genKstTrkmPx',           genKstTrkmPx,                           'genKstTrkmPx/D')
# ntuple.Branch('genKstTrkmPy',           genKstTrkmPy,                           'genKstTrkmPy/D')
# ntuple.Branch('genKstTrkmPz',           genKstTrkmPz,                           'genKstTrkmPz/D')
ntuple.Branch('genKstTrkmID',           genKstTrkmID,                           'genKstTrkmID/D')
ntuple.Branch('genKstTrkpMother',       genKstTrkpMother,                       'genKstTrkpMother/D')
# ntuple.Branch('genKstTrkpPx',           genKstTrkpPx,                           'genKstTrkpPx/D')
# ntuple.Branch('genKstTrkpPy',           genKstTrkpPy,                           'genKstTrkpPy/D')
# ntuple.Branch('genKstTrkpPz',           genKstTrkpPz,                           'genKstTrkpPz/D')
ntuple.Branch('genKstTrkpID',           genKstTrkpID,                           'genKstTrkpID/D')

ntuple.Branch('genbPt',        genbPt          ,               'genbPt/D')
ntuple.Branch('genkstPt',      genkstPt        ,               'genkstPt/D')
# ntuple.Branch('genmumuPt',     genmumuPt       ,               'genmumuPt/D')
ntuple.Branch('genmumPt',      genmumPt        ,               'genmumPt/D')
ntuple.Branch('genmupPt',      genmupPt        ,               'genmupPt/D')
ntuple.Branch('genkstTrkmPt',  genkstTrkmPt    ,               'genkstTrkmPt/D')
ntuple.Branch('genkstTrkpPt',  genkstTrkpPt    ,               'genkstTrkpPt/D')
ntuple.Branch('genbPhi',       genbPhi         ,               'genbPhi/D')
ntuple.Branch('genkstPhi',     genkstPhi       ,               'genkstPhi/D')
# ntuple.Branch('genmumuPhi',    genmumuPhi      ,               'genmumuPhi/D')
ntuple.Branch('genmumPhi',     genmumPhi       ,               'genmumPhi/D')
ntuple.Branch('genmupPhi',     genmupPhi       ,               'genmupPhi/D')
ntuple.Branch('genkstTrkmPhi', genkstTrkmPhi   ,               'genkstTrkmPhi/D')
ntuple.Branch('genkstTrkpPhi', genkstTrkpPhi   ,               'genkstTrkpPhi/D')
ntuple.Branch('genbEta',       genbEta         ,               'genbEta/D')
ntuple.Branch('genkstEta',     genkstEta       ,               'genkstEta/D')
# ntuple.Branch('genmumuEta',    genmumuEta    ,               'genmumuEta/D')
ntuple.Branch('genmumEta',     genmumEta       ,               'genmumEta/D')
ntuple.Branch('genmupEta',     genmupEta       ,               'genmupEta/D')
ntuple.Branch('genkstTrkmEta', genkstTrkmEta   ,               'genkstTrkmEta/D')
ntuple.Branch('genkstTrkpEta', genkstTrkpEta   ,               'genkstTrkpEta/D')

ntuple.Branch('cos_theta_l' ,  cos_theta_l  ,               'cos_theta_l/D')
ntuple.Branch('cos_theta_k' ,  cos_theta_k  ,               'cos_theta_k/D')
ntuple.Branch('phi_kst_mumu',  phi_kst_mumu ,               'phi_kst_mumu/D')

ntuple.Branch('genQ2',         genQ2           ,               'genQ2/D')
ntuple.Branch('genQ',          genQ            ,               'genQ/D')


numEvents = tree_lmnr.GetEntries()
print 'total number of events in tree:', numEvents


progressbarWidth = 40
sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
sys.stdout.flush()                          # this forces to print the stdout buffer
sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['


tree_lmnr.SetBranchAddress('eventN',eventN)
for i, ev in enumerate(tree_lmnr):

    if i%int(numEvents/(progressbarWidth-1))==0:
        sys.stdout.write('+')
        sys.stdout.flush()

#     if i > 10000: 
#         break

    for var in genbr:
	var[0] = -99.

    ## per event quantities
    runN[0]                        = ev.runN

    ## per event GEN quantities
    genSignal[0]                   = ev.genSignal
    genSignHasFSR[0]               = ev.genSignHasFSR
    genSignKstHasFSR[0]            = ev.genSignKstHasFSR
    genSignPsiHasFSR[0]            = ev.genSignPsiHasFSR
#     genPriVtxX[0]                  = ev.genPriVtxX
#     genPriVtxY[0]                  = ev.genPriVtxY
#     genPriVtxZ[0]                  = ev.genPriVtxZ
    genB0Mass[0]                   = ev.genB0Mass
#     genB0Px[0]                     = ev.genB0Px
#     genB0Py[0]                     = ev.genB0Py
#     genB0Pz[0]                     = ev.genB0Pz
#     genB0VtxX[0]                   = ev.genB0VtxX
#     genB0VtxY[0]                   = ev.genB0VtxY
#     genB0VtxZ[0]                   = ev.genB0VtxZ
    genKstMass[0]                  = ev.genKstMass
#     genKstPx[0]                    = ev.genKstPx
#     genKstPy[0]                    = ev.genKstPy
#     genKstPz[0]                    = ev.genKstPz
#     genKstVtxX[0]                  = ev.genKstVtxX
#     genKstVtxY[0]                  = ev.genKstVtxY
#     genKstVtxZ[0]                  = ev.genKstVtxZ
    genPsiMass[0]                  = ev.genPsiMass
#     genPsiVtxX[0]                  = ev.genPsiVtxX
#     genPsiVtxY[0]                  = ev.genPsiVtxY
#     genPsiVtxZ[0]                  = ev.genPsiVtxZ
    genMumMother[0]                = ev.genMumMother
#     genMumPx[0]                    = ev.genMumPx
#     genMumPy[0]                    = ev.genMumPy
#     genMumPz[0]                    = ev.genMumPz
    genMupMother[0]                = ev.genMupMother
#     genMupPx[0]                    = ev.genMupPx
#     genMupPy[0]                    = ev.genMupPy
#     genMupPz[0]                    = ev.genMupPz
    genKstTrkmMother[0]            = ev.genKstTrkmMother
    genKstTrkmID[0]                = ev.genKstTrkmID
#     genKstTrkmPx[0]                = ev.genKstTrkmPx
#     genKstTrkmPy[0]                = ev.genKstTrkmPy
#     genKstTrkmPz[0]                = ev.genKstTrkmPz
    genKstTrkpMother[0]            = ev.genKstTrkpMother
    genKstTrkpID[0]                = ev.genKstTrkpID
#     genKstTrkpPx[0]                = ev.genKstTrkpPx
#     genKstTrkpPy[0]                = ev.genKstTrkpPy
#     genKstTrkpPz[0]                = ev.genKstTrkpPz

    genbPt[0]                      = computePt ( ev.genB0Px,  ev.genB0Py ) 
    genkstPt[0]                    = computePt ( ev.genKstPx, ev.genKstPy ) 
#     genmumuPt[0]                   = computePt( ev.mumuPx   [icand], ev.mumuPy    [icand] ) 
    genmumPt[0]                    = computePt ( ev.genMumPx, ev.genMumPy) 
    genmupPt[0]                    = computePt ( ev.genMupPx, ev.genMupPy) 
    genkstTrkmPt[0]                = computePt ( ev.genKstTrkmPx, ev.genKstTrkmPy ) 
    genkstTrkpPt[0]                = computePt ( ev.genKstTrkpPx, ev.genKstTrkpPy ) 
    genbPhi[0]                     = computePhi( ev.genB0Px,  ev.genB0Py )
    genkstPhi[0]                   = computePhi( ev.genKstPx, ev.genKstPy )
#     genmumuPhi[0]                  = computePhi( ev.mumuPx   [icand], ev.mumuPy    [icand] )
    genmumPhi[0]                   = computePhi( ev.genMumPx, ev.genMumPy)
    genmupPhi[0]                   = computePhi( ev.genMupPx, ev.genMupPy)
    genkstTrkmPhi[0]               = computePhi( ev.genKstTrkmPx, ev.genKstTrkmPy )
    genkstTrkpPhi[0]               = computePhi( ev.genKstTrkpPx, ev.genKstTrkpPy )
    genbEta[0]                     = computeEta( ev.genB0Px      , ev.genB0Py       , ev.genB0Pz        )
    genkstEta[0]                   = computeEta( ev.genKstPx    , ev.genKstPy     , ev.genKstPz      )
#     genmumuEta[0]                  = computeEta( ev.mumuPx   , ev.mumuPy    , ev.mumuPz     )
    genmumEta[0]                   = computeEta( ev.genMumPx    , ev.genMumPy     , ev.genMumPz      )
    genmupEta[0]                   = computeEta( ev.genMupPx    , ev.genMupPy     , ev.genMupPz      )
    genkstTrkmEta[0]               = computeEta( ev.genKstTrkmPx, ev.genKstTrkmPy , ev.genKstTrkmPz  )
    genkstTrkpEta[0]               = computeEta( ev.genKstTrkpPx, ev.genKstTrkpPy , ev.genKstTrkpPz  )
    genQ[0]                        = computeInvMass(ev.genMumPx, ev.genMumPy, ev.genMumPz, muonMass,
                                                    ev.genMupPx, ev.genMupPy, ev.genMupPz, muonMass )
    genQ2[0]                       = genQ[0]*genQ[0]
                                                                                                            
                                                    
    cos_theta_l[0], cos_theta_k[0], phi_kst_mumu[0] = addGenVars(
    	ev.genSignal,
    	ev.genMumPx,      ev.genMumPy,     ev.genMumPz,  
    	ev.genMupPx,      ev.genMupPy,     ev.genMupPz,  
    	ev.genKstTrkmPx,  ev.genKstTrkmPy, ev.genKstTrkmPz,
    	ev.genKstTrkpPx,  ev.genKstTrkpPy, ev.genKstTrkpPz,
    	ev.genKstPx,      ev.genKstPy,     ev.genKstPz,   ev.genKstMass,
    	ev.genB0Px,       ev.genB0Py,      ev.genB0Pz,    ev.genB0Mass
    )
                                                    



    ntuple. Fill()


sys.stdout.write('\n')

file_out.cd()
ntuple.Write()
file_out.Close()
