import ROOT
import sys
import math
from DataFormats.FWLite import Events
from collections import OrderedDict
from array import array
import numpy
from ROOT import TLorentzVector
import cProfile, pstats, io

'''
'''

pionmass_ = 0.139570

b0_mass  = 5.27963
b0_width = 0.03601 
kaonMass = 0.493677
muonMass = 0.10565837

kstMass_     = 0.896
kstSigma_    = 0.05

type = 'LMNR'
tree_lmnr = ROOT.TChain('B0KstMuMu/B0KstMuMuNTuple')

## Jpsi MC
if type == 'JPsi':
   tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_B0KstarJpsi_addIso_sub2/180116_054002/0002/B0ToKstMuMu_*.root')

## low mass non resonant
if type == 'LMNR':
 	tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/phaseII/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_PhaseII_KStarMuMu_firstRound_may30/180530_122915/0000/B0ToKstMuMu_*.root')


file_out = ROOT.TFile('ntuples/flat_GEN_LMNR.root', 'recreate')
# file_out = ROOT.TFile('ntuples/mc/out_flat_2016MC_allGEN_sub3_p1.root', 'recreate')
ntuple   = ROOT.TTree( 'ntuple', 'ntuple' )



def computeEta(px,py,pz):
    P = math.sqrt(px*px + py*py + pz*pz)
    return 0.5*math.log((P + pz) / (P - pz))
    
def computePhi(px,py):
    phi = math.atan(py / px)
    if (px < 0 and py < 0):  phi = phi - math.pi;
    if (px < 0 and py > 0):  phi = phi + math.pi;
    return phi

def computePt(px,py):
     return math.sqrt(px*px+py*py) 
 
def computeInvMass(px1,py1,pz1,m1,px2,py2,pz2,m2):
     e1 = math.sqrt(px1**2 + py1**2 + pz1**2 + m1**2)
     e2 = math.sqrt(px2**2 + py2**2 + pz2**2 + m2**2)
     return math.sqrt( (e1+e2)**2 - ((px1+px2)**2 + (py1+py2)**2 + (pz1+pz2)**2 ) )
 
def computeCosine (Vx, Vy, Vz,
			       Wx, Wy, Wz):

  Vnorm = math.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
  Wnorm = math.sqrt(Wx*Wx + Wy*Wy + Wz*Wz)
  VdotW = Vx*Wx + Vy*Wy + Vz*Wz
  
  if Vnorm > 0. and Wnorm > 0.:
    cosAlpha = VdotW / (Vnorm * Wnorm)
  else:
    cosAlpha = -99.
    
  return cosAlpha


mumu_lv         = TLorentzVector()
kst_lv          = TLorentzVector()
b0_lv           = TLorentzVector()
# mum_lv          = TLorentzVector()
# mup_lv          = TLorentzVector()
# tkp_lv          = TLorentzVector()
# tkm_lv          = TLorentzVector()
mu_lv_boosted   = TLorentzVector()
kaon_lv_boosted = TLorentzVector()
b0_lv_boosted   = TLorentzVector()
mum_lv_boosted  = TLorentzVector()
mup_lv_boosted  = TLorentzVector()
tkp_lv_boosted  = TLorentzVector()
tkm_lv_boosted  = TLorentzVector()

def addGenVars(genSignal,
           mumPx,      mumPy,     mumPz,  
           mupPx,      mupPy,     mupPz,  
           kstTrkmPx,  kstTrkmPy, kstTrkmPz,
           kstTrkpPx,  kstTrkpPy, kstTrkpPz,
           kstPx,      kstPy,     kstPz,   kstMass,
           bPx,        bPy,       bPz,     bMass
          ):

  tagB0 = 1 if genSignal==1 else 0        
  
  kst_e  = math.sqrt( kstPx**2 + kstPy**2 + kstPz**2 + kstMass**2 )
  kst_lv.SetPxPyPzE(kstPx, kstPy, kstPz, kst_e)

  b0_e  = math.sqrt( bPx**2 + bPy**2 + bPz**2 + bMass**2 )
  b0_lv.SetPxPyPzE(bPx, bPy, bPz, b0_e)
  
  mum_e  = math.sqrt( mumPx**2 + mumPy**2 + mumPz**2 + muonMass**2 )
  mup_e  = math.sqrt( mupPx**2 + mupPy**2 + mupPz**2 + muonMass**2 )
#   mum_lv.SetPxPyPzE(mumPx, mumPy, mumPz, mum_e)
#   mup_lv.SetPxPyPzE(mupPx, mupPy, mupPz, mup_e)
#   mumu_lv = mup_lv + mum_lv

  mumu_lv.SetPxPyPzE(mumPx + mupPx, mumPy + mupPy, mumPz + mupPz, mum_e + mup_e)


  tkp_e  = math.sqrt( kstTrkpPx**2 + kstTrkpPy**2 + kstTrkpPz**2 + tagB0*kaonMass**2  + (1-tagB0)*pionmass_**2 )
  tkm_e  = math.sqrt( kstTrkmPx**2 + kstTrkmPy**2 + kstTrkmPz**2 + tagB0*pionmass_**2 + (1-tagB0)*kaonMass**2  )
  
  cos_theta_k     = -99.
  cos_theta_l     = -99.
  phiKstMuMuPlane = -99.

  if tagB0:
    mu_lv_boosted  .SetPxPyPzE(mupPx, mupPy, mupPz, mup_e)
    kaon_lv_boosted.SetPxPyPzE(kstTrkpPx,kstTrkpPy,kstTrkpPz,tkp_e)
  
  else:
    mu_lv_boosted  .SetPxPyPzE(mumPx, mumPy, mumPz, mum_e)
    kaon_lv_boosted.SetPxPyPzE(kstTrkmPx, kstTrkmPy, kstTrkmPz, tkm_e)
    
  ## calculate cos theta_l
  boostMuMu = mumu_lv.BoostVector()
  mu_lv_boosted.Boost(-boostMuMu)
  b0_lv_boosted  .SetPxPyPzE(b0_lv.Px(),b0_lv.Py(),b0_lv.Pz(),b0_lv.E())
  b0_lv_boosted.Boost(-boostMuMu)
  
  cos_theta_l = computeCosine(-b0_lv_boosted.Px(),
                              -b0_lv_boosted.Py(),
                              -b0_lv_boosted.Pz(),
                               mu_lv_boosted.Px(),
                               mu_lv_boosted.Py(),
                               mu_lv_boosted.Pz()
                             )

  ## calculate cos theta_k
  boostKst = kst_lv.BoostVector()
  kaon_lv_boosted.Boost(-boostKst)
  b0_lv_boosted  .SetPxPyPzE(b0_lv.Px(),b0_lv.Py(),b0_lv.Pz(),b0_lv.E())
  b0_lv_boosted.Boost(-boostKst)

  cos_theta_k = computeCosine(-b0_lv_boosted.Px(),
                              -b0_lv_boosted.Py(),
                              -b0_lv_boosted.Pz(),
                               kaon_lv_boosted.Px(),
                               kaon_lv_boosted.Py(),
                               kaon_lv_boosted.Pz()
                             )

  ## calculate angle between planes
  boostB0 = b0_lv.BoostVector()
  mum_lv_boosted  .SetPxPyPzE(mumPx, mumPy, mumPz, mum_e)
  mum_lv_boosted.Boost(-boostB0)

  mup_lv_boosted  .SetPxPyPzE(mupPx, mupPy, mupPz, mup_e)
  mup_lv_boosted.Boost(-boostB0)

  tkp_lv_boosted  .SetPxPyPzE(kstTrkpPx, kstTrkpPy, kstTrkpPz, tkp_e)
  tkp_lv_boosted.Boost(-boostB0)

  tkm_lv_boosted  .SetPxPyPzE(kstTrkmPx, kstTrkmPy, kstTrkmPz, tkm_e)
  tkm_lv_boosted.Boost(-boostB0)
  
  MuMuPlane = mup_lv_boosted.Vect().Cross(mum_lv_boosted.Vect())
  KstPlane  = tkp_lv_boosted.Vect().Cross(tkm_lv_boosted.Vect())

  phiKstMuMuPlane = -99
  if MuMuPlane.Cross(KstPlane).Dot(-b0_lv.Vect()) > 0:
    phiKstMuMuPlane = MuMuPlane.Angle(KstPlane)
  else:
    phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane)

  return cos_theta_l, cos_theta_k, phiKstMuMuPlane



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
runN                             = array('d', [-99.]);  genbr.append(runN)
eventN                           = array('d', [-99.]);  genbr.append(eventN)
genSignal                        = array('d', [-99.]);  genbr.append(genSignal)
genSignHasFSR                    = array('d', [-99.]);  genbr.append(genSignHasFSR)
genSignKstHasFSR                 = array('d', [-99.]);  genbr.append(genSignKstHasFSR)
genSignPsiHasFSR                 = array('d', [-99.]);  genbr.append(genSignPsiHasFSR)
genPriVtxX                       = array('d', [-99.]);  genbr.append(genPriVtxX)
genPriVtxY                       = array('d', [-99.]);  genbr.append(genPriVtxY)
genPriVtxZ                       = array('d', [-99.]);  genbr.append(genPriVtxZ)
genB0Mass                        = array('d', [-99.]);  genbr.append(genB0Mass)
genB0Px                          = array('d', [-99.]);  genbr.append(genB0Px)
genB0Py                          = array('d', [-99.]);  genbr.append(genB0Py)
genB0Pz                          = array('d', [-99.]);  genbr.append(genB0Pz)
genB0VtxX                        = array('d', [-99.]);  genbr.append(genB0VtxX)
genB0VtxY                        = array('d', [-99.]);  genbr.append(genB0VtxY)
genB0VtxZ                        = array('d', [-99.]);  genbr.append(genB0VtxZ)
genKstMass                       = array('d', [-99.]);  genbr.append(genKstMass)
genKstPx                         = array('d', [-99.]);  genbr.append(genKstPx)
genKstPy                         = array('d', [-99.]);  genbr.append(genKstPy)
genKstPz                         = array('d', [-99.]);  genbr.append(genKstPz)
genKstVtxX                       = array('d', [-99.]);  genbr.append(genKstVtxX)
genKstVtxY                       = array('d', [-99.]);  genbr.append(genKstVtxY)
genKstVtxZ                       = array('d', [-99.]);  genbr.append(genKstVtxZ)
genPsiMass                       = array('d', [-99.]);  genbr.append(genPsiMass)
genPsiVtxX                       = array('d', [-99.]);  genbr.append(genPsiVtxX)
genPsiVtxY                       = array('d', [-99.]);  genbr.append(genPsiVtxY)
genPsiVtxZ                       = array('d', [-99.]);  genbr.append(genPsiVtxZ)
genMumMother                     = array('d', [-99.]);  genbr.append(genMumMother)
genMumPx                         = array('d', [-99.]);  genbr.append(genMumPx)
genMumPy                         = array('d', [-99.]);  genbr.append(genMumPy)
genMumPz                         = array('d', [-99.]);  genbr.append(genMumPz)
genMupMother                     = array('d', [-99.]);  genbr.append(genMupMother)
genMupPx                         = array('d', [-99.]);  genbr.append(genMupPx)
genMupPy                         = array('d', [-99.]);  genbr.append(genMupPy)
genMupPz                         = array('d', [-99.]);  genbr.append(genMupPz)
genKstTrkmMother                 = array('d', [-99.]);  genbr.append(genKstTrkmMother)
genKstTrkmID                     = array('d', [-99.]);  genbr.append(genKstTrkmID)
genKstTrkmPx                     = array('d', [-99.]);  genbr.append(genKstTrkmPx)
genKstTrkmPy                     = array('d', [-99.]);  genbr.append(genKstTrkmPy)
genKstTrkmPz                     = array('d', [-99.]);  genbr.append(genKstTrkmPz)
genKstTrkpMother                 = array('d', [-99.]);  genbr.append(genKstTrkpMother)
genKstTrkpID                     = array('d', [-99.]);  genbr.append(genKstTrkpID)
genKstTrkpPx                     = array('d', [-99.]);  genbr.append(genKstTrkpPx)
genKstTrkpPy                     = array('d', [-99.]);  genbr.append(genKstTrkpPy)
genKstTrkpPz                     = array('d', [-99.]);  genbr.append(genKstTrkpPz)

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
genq2                               = array('d',[-99.]); genbr.append( genq2           )   

gen_cos_theta_l                     = array('d',[-99.]); genbr.append( gen_cos_theta_l )   
gen_cos_theta_k                     = array('d',[-99.]); genbr.append( gen_cos_theta_k )   
gen_phi_kst_mumu                    = array('d',[-99.]); genbr.append( gen_phi_kst_mumu)   


'''
for i in branches:
    print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/D')".format(NAME=i).expandtabs(40)
'''
ntuple.Branch('runN',                   runN,                                   'runN/D')
ntuple.Branch('eventN',                 eventN,                                 'eventN/D')
ntuple.Branch('genSignal',              genSignal,                              'genSignal/D')
ntuple.Branch('genSignHasFSR',          genSignHasFSR,                          'genSignHasFSR/D')
ntuple.Branch('genSignKstHasFSR',       genSignKstHasFSR,                       'genSignKstHasFSR/D')
ntuple.Branch('genSignPsiHasFSR',       genSignPsiHasFSR,                       'genSignPsiHasFSR/D')
ntuple.Branch('genPriVtxX',             genPriVtxX,                             'genPriVtxX/D')
ntuple.Branch('genPriVtxY',             genPriVtxY,                             'genPriVtxY/D')
ntuple.Branch('genPriVtxZ',             genPriVtxZ,                             'genPriVtxZ/D')
ntuple.Branch('genB0Mass',              genB0Mass,                              'genB0Mass/D')
ntuple.Branch('genB0Px',                genB0Px,                                'genB0Px/D')
ntuple.Branch('genB0Py',                genB0Py,                                'genB0Py/D')
ntuple.Branch('genB0Pz',                genB0Pz,                                'genB0Pz/D')
ntuple.Branch('genB0VtxX',              genB0VtxX,                              'genB0VtxX/D')
ntuple.Branch('genB0VtxY',              genB0VtxY,                              'genB0VtxY/D')
ntuple.Branch('genB0VtxZ',              genB0VtxZ,                              'genB0VtxZ/D')
ntuple.Branch('genKstMass',             genKstMass,                             'genKstMass/D')
ntuple.Branch('genKstPx',               genKstPx,                               'genKstPx/D')
ntuple.Branch('genKstPy',               genKstPy,                               'genKstPy/D')
ntuple.Branch('genKstPz',               genKstPz,                               'genKstPz/D')
ntuple.Branch('genKstVtxY',             genKstVtxY,                             'genKstVtxY/D')
ntuple.Branch('genPsiMass',             genPsiMass,                             'genPsiMass/D')
ntuple.Branch('genPsiVtxX',             genPsiVtxX,                             'genPsiVtxX/D')
ntuple.Branch('genPsiVtxY',             genPsiVtxY,                             'genPsiVtxY/D')
ntuple.Branch('genPsiVtxZ',             genPsiVtxZ,                             'genPsiVtxZ/D')
ntuple.Branch('genMumMother',           genMumMother,                           'genMumMother/D')
ntuple.Branch('genMumPx',               genMumPx,                               'genMumPx/D')
ntuple.Branch('genMumPy',               genMumPy,                               'genMumPy/D')
ntuple.Branch('genMumPz',               genMumPz,                               'genMumPz/D')
ntuple.Branch('genMupMother',           genMupMother,                           'genMupMother/D')
ntuple.Branch('genMupPx',               genMupPx,                               'genMupPx/D')
ntuple.Branch('genMupPy',               genMupPy,                               'genMupPy/D')
ntuple.Branch('genMupPz',               genMupPz,                               'genMupPz/D')
ntuple.Branch('genKstTrkmMother',       genKstTrkmMother,                       'genKstTrkmMother/D')
ntuple.Branch('genKstTrkmPx',           genKstTrkmPx,                           'genKstTrkmPx/D')
ntuple.Branch('genKstTrkmPy',           genKstTrkmPy,                           'genKstTrkmPy/D')
ntuple.Branch('genKstTrkmPz',           genKstTrkmPz,                           'genKstTrkmPz/D')
ntuple.Branch('genKstTrkmID',           genKstTrkmID,                           'genKstTrkmID/D')
ntuple.Branch('genKstTrkpMother',       genKstTrkpMother,                       'genKstTrkpMother/D')
ntuple.Branch('genKstTrkpPx',           genKstTrkpPx,                           'genKstTrkpPx/D')
ntuple.Branch('genKstTrkpPy',           genKstTrkpPy,                           'genKstTrkpPy/D')
ntuple.Branch('genKstTrkpPz',           genKstTrkpPz,                           'genKstTrkpPz/D')
ntuple.Branch('genKstTrkpID',           genKstTrkpID,                           'genKstTrkpID/D')
ntuple.Branch('genq2',                  genq2,                                  'genq2/D')

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

ntuple.Branch('gen_cos_theta_l' , gen_cos_theta_l  ,               'gen_cos_theta_l/D')
ntuple.Branch('gen_cos_theta_k' , gen_cos_theta_k  ,               'gen_cos_theta_k/D')
ntuple.Branch('gen_phi_kst_mumu', gen_phi_kst_mumu ,               'gen_phi_kst_mumu/D')


numEvents = tree_lmnr.GetEntries()
print 'total number of events in tree:', numEvents


progressbarWidth = 40
sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
sys.stdout.flush()                          # this forces to print the stdout buffer
sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['

pr = cProfile.Profile()
pr.enable()



for i, ev in enumerate(tree_lmnr):

    if i%int(numEvents/(progressbarWidth-1))==0:
        sys.stdout.write('+')
        sys.stdout.flush()

#     if i > 50000: 
#         break

	for var in genbr:
		var[0] = -99.

    ## per event quantities
    runN[0]                        = ev.runN
    eventN[0]                      = ev.eventN

    ## per event GEN quantities
    genSignal[0]                   = ev.genSignal
    genSignHasFSR[0]               = ev.genSignHasFSR
    genSignKstHasFSR[0]            = ev.genSignKstHasFSR
    genSignPsiHasFSR[0]            = ev.genSignPsiHasFSR
    genPriVtxX[0]                  = ev.genPriVtxX
    genPriVtxY[0]                  = ev.genPriVtxY
    genPriVtxZ[0]                  = ev.genPriVtxZ
    genB0Mass[0]                   = ev.genB0Mass
    genB0Px[0]                     = ev.genB0Px
    genB0Py[0]                     = ev.genB0Py
    genB0Pz[0]                     = ev.genB0Pz
    genB0VtxX[0]                   = ev.genB0VtxX
    genB0VtxY[0]                   = ev.genB0VtxY
    genB0VtxZ[0]                   = ev.genB0VtxZ
    genKstMass[0]                  = ev.genKstMass
    genKstPx[0]                    = ev.genKstPx
    genKstPy[0]                    = ev.genKstPy
    genKstPz[0]                    = ev.genKstPz
    genKstVtxX[0]                  = ev.genKstVtxX
    genKstVtxY[0]                  = ev.genKstVtxY
    genKstVtxZ[0]                  = ev.genKstVtxZ
    genPsiMass[0]                  = ev.genPsiMass
    genPsiVtxX[0]                  = ev.genPsiVtxX
    genPsiVtxY[0]                  = ev.genPsiVtxY
    genPsiVtxZ[0]                  = ev.genPsiVtxZ
    genMumMother[0]                = ev.genMumMother
    genMumPx[0]                    = ev.genMumPx
    genMumPy[0]                    = ev.genMumPy
    genMumPz[0]                    = ev.genMumPz
    genMupMother[0]                = ev.genMupMother
    genMupPx[0]                    = ev.genMupPx
    genMupPy[0]                    = ev.genMupPy
    genMupPz[0]                    = ev.genMupPz
    genKstTrkmMother[0]            = ev.genKstTrkmMother
    genKstTrkmID[0]                = ev.genKstTrkmID
    genKstTrkmPx[0]                = ev.genKstTrkmPx
    genKstTrkmPy[0]                = ev.genKstTrkmPy
    genKstTrkmPz[0]                = ev.genKstTrkmPz
    genKstTrkpMother[0]            = ev.genKstTrkpMother
    genKstTrkpID[0]                = ev.genKstTrkpID
    genKstTrkpPx[0]                = ev.genKstTrkpPx
    genKstTrkpPy[0]                = ev.genKstTrkpPy
    genKstTrkpPz[0]                = ev.genKstTrkpPz

    genbPt[0]                      = computePt ( ev.genB0Px,  ev.genB0Py ) 
    genkstPt[0]                    = computePt ( ev.genKstPx, ev.genKstPy ) 
#     genmumuPt[0]                   = computePt( ev.mumuPx   [icand], ev.mumuPy    [icand] ) 
    genmumPt[0]                    = computePt ( ev.genMumPx, ev.genMumPy) 
    genmupPt[0]                    = computePt ( ev.genMupPx, ev.genMupPx) 
    genkstTrkmPt[0]                = computePt ( ev.genKstTrkmPx, ev.genKstTrkmPy ) 
    genkstTrkpPt[0]                = computePt ( ev.genKstTrkpPx, ev.genKstTrkpPy ) 
    genbPhi[0]                     = computePhi( ev.genB0Px,  ev.genB0Py )
    genkstPhi[0]                   = computePhi( ev.genKstPx, ev.genKstPy )
#     genmumuPhi[0]                  = computePhi( ev.mumuPx   [icand], ev.mumuPy    [icand] )
    genmumPhi[0]                   = computePhi( ev.genMumPx, ev.genMumPy)
    genmupPhi[0]                   = computePhi( ev.genMupPx, ev.genMupPx)
    genkstTrkmPhi[0]               = computePhi( ev.genKstTrkmPx, ev.genKstTrkmPy )
    genkstTrkpPhi[0]               = computePhi( ev.genKstTrkpPx, ev.genKstTrkpPy )
    genbEta[0]                     = computeEta( ev.genB0Px      , ev.genB0Py       , ev.genB0Pz        )
    genkstEta[0]                   = computeEta( ev.genKstPx    , ev.genKstPy     , ev.genKstPz      )
#     genmumuEta[0]                  = computeEta( ev.mumuPx   , ev.mumuPy    , ev.mumuPz     )
    genmumEta[0]                   = computeEta( ev.genMumPx    , ev.genMumPy     , ev.genMumPz      )
    genmupEta[0]                   = computeEta( ev.genMupPx    , ev.genMupPy     , ev.genMupPz      )
    genkstTrkmEta[0]               = computeEta( ev.genKstTrkmPx, ev.genKstTrkmPy , ev.genKstTrkmPz  )
    genkstTrkpEta[0]               = computeEta( ev.genKstTrkpPx, ev.genKstTrkpPy , ev.genKstTrkpPz  )
    genq2[0]                       = computeInvMass(ev.genMumPx, ev.genMumPy, ev.genMumPz, muonMass,
                                                    ev.genMupPx, ev.genMupPy, ev.genMupPz, muonMass )
                                                    
                                                    
    gen_cos_theta_l[0], gen_cos_theta_k[0], gen_phi_kst_mumu[0] = addGenVars(
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

pr.disable()
s = open('x.prof', 'a')
sortby = 'cumulative'
pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby).print_stats()
s.close()










'''
branches = [
            'runN',
            'eventN',
            'genSignal',        
            'genSignHasFSR',    
            'genSignKstHasFSR', 
            'genSignPsiHasFSR', 
            'genPriVtxX',       
            'genPriVtxY',       
            'genPriVtxZ',       
            'genB0Mass',        
            'genB0Px',          
            'genB0Py',          
            'genB0Pz',          
            'genB0VtxX',        
            'genB0VtxY',        
            'genB0VtxZ',        
            'genKstMass',       
            'genKstPx',         
            'genKstPy',         
            'genKstPz',         
            'genKstVtxX',       
            'genKstVtxY',       
            'genKstVtxZ',       
            'genPsiMass',       
            'genPsiVtxX',       
            'genPsiVtxY',       
            'genPsiVtxZ',       
            'genMumMother',     
            'genMumPx',         
            'genMumPy',         
            'genMumPz',         
            'genMupMother',     
            'genMupPx',         
            'genMupPy',         
            'genMupPz',         
            'genKstTrkmMother', 
            'genKstTrkmID',     
            'genKstTrkmPx',     
            'genKstTrkmPy',     
            'genKstTrkmPz',     
            'genKstTrkpMother', 
            'genKstTrkpID',     
            'genKstTrkpPx',     
            'genKstTrkpPy',     
            'genKstTrkpPz'
]
'''

