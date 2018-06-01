import ROOT
import sys
import math
from DataFormats.FWLite import Events
from collections import OrderedDict
from array import array
import numpy

'''
STILL TODO: hadd info about which trigger the muons and tracks are matched to
(for current 2016 ntuples not needed because 1 trigger per crab job)
'''


b0_mass  = 5.27963
b0_width = 0.03601 
kaonMass = 0.493677

kstMass_     = 0.896
kstSigma_    = 0.05
KstMassShape = ROOT.TF1('KstMassShape', '2*sqrt(2)*[0]*[1]* sqrt([0]*[0] * ([0]*[0] + [1]*[1])) / (TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0] * ([0]*[0] + [1]*[1])))) / ((x*x - [0]*[0]) * (x*x - [0]*[0]) + [0]*[0]*[1]*[1])', 0.0,kstMass_*2.);
KstMassShape.SetParameter(0,kstMass_)
KstMassShape.SetParameter(1,kstSigma_)
KstMassShape.SetParNames("Mean","Width")

###############################################
## Parameters for the flattening:

paths = [ 
          'HLT_DoubleMu4_LowMassNonResonantTrk_Displaced',
          'HLT_DoubleMu4_JpsiTrk_Displaced',
          'HLT_DoubleMu4_PsiPrimeTrk_Displaced',
        ]  


type       = 'LMNR'
skim       = False
skimSoftMu = False
###############################################


categories = [
           'GlobalMuon ',   
           'TrackerMuon ',   
           'StandAloneMuon',   
           'TrackerMuonArbitrated ',   
           'GlobalMuonPromptTight ',   
           'TMLastStationTight ',   
           'TMLastStationLoose ',   
           'TM2DCompatibilityTight ',   
           'TM2DCompatibilityLoose ',   
           'TMOneStationTight ',   
           'TMOneStationLoose ',   
           'TMLastStationAngTight ',   
           'TMLastStationAngLoose ',   
           'TMOneStationAngTight ',   
           'TMOneStationAngLoose '
          ]   

tree_lmnr = ROOT.TChain('B0KstMuMu/B0KstMuMuNTuple')

## Jpsi MC
if type == 'JPsi':
  tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_B0KstarJpsi_addIso_sub1/180115_204406/0000/B0ToKstMuMu_*.root')

## low mass non resonant
if type == 'LMNR':
 	tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/phaseII/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_PhaseII_KStarMuMu_firstRound_may30/180530_122915/0000/B0ToKstMuMu_*.root')

## Mauro 2012, GEN
# tree_lmnr.Add('/gwteras/cms/store/user/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/GENcands/B0ToKstMuMu_GEN_NoFilter_MC_NTuples.root')


file_out = ROOT.TFile('ntuples/flat_ntuple_PhaseII_RECO_%s.root'%type, 'recreate')
ntuple   = ROOT.TTree( 'ntuple', 'ntuple' )

isMC = False
if 'MC' in tree_lmnr.GetFile().GetName() or 'pythia' in tree_lmnr.GetFile().GetName():
    isMC = True



def FlavorTagger(thekstMass,thekstBarMass):
    ''' choose B0 or B0bar based on the K+/pi- K-/pi+ invariant mass'''
    distKst       = abs(thekstMass    - kstMass_)
    distKstBar    = abs(thekstBarMass - kstMass_)
    return 1 if distKst < distKstBar else 0


def muonCategory(mucat):
    '''unpack muon category from original ntuples'''

    muCategoryDict = {}
    for cat in categories:
        muCategoryDict[cat.strip()] = 1 if cat in mucat else 0
    return muCategoryDict       

def muonIso(muiso, num=20):
    '''calculate # iso tracks below threshold'''

    counts = numpy.zeros(num)   
    if len(muiso) > 0:
        for i in range(num):
            counts[i] = numpy.size(numpy.where(numpy.array(muiso) < 0.1*(i+1)))

    toreturn = ROOT.std.vector('float')() 
    for i in counts:  toreturn.push_back(i)
    return toreturn


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
 




print '@@@@@@@@@@@     flattening ntuple     @@@@@@@@@@@'
print 'input file: ', tree_lmnr.GetFile()
print 'Dataset:', isMC*'MC', (not isMC)*'Data'

try:
  tree_lmnr.GetEntries()
except:
  print 'ntuple tree not found'
  exit()


newbr =  []

'''
for i in branches:
    print "{NAME}\t = array('f', [-99.]);  newbr.append({NAME})".format(NAME=i).expandtabs(40)
'''
runN                             = array('f', [-99.]);  newbr.append(runN)
eventN                           = array('f', [-99.]);  newbr.append(eventN)
recoVtxN                         = array('f', [-99.]);  newbr.append(recoVtxN)
evWeight                         = array('f', [-99.]);  newbr.append(evWeight)
evWeightE2                       = array('f', [-99.]);  newbr.append(evWeightE2)
trueNumInteractionsMC            = array('f', [-99.]);  newbr.append(trueNumInteractionsMC)
trig                             = array('f', [-99.]);  newbr.append(trig)
bsX                              = array('f', [-99.]);  newbr.append(bsX)
bsY                              = array('f', [-99.]);  newbr.append(bsY)
bMass                            = array('f', [-99.]);  newbr.append(bMass)
bMassE                           = array('f', [-99.]);  newbr.append(bMassE)
bBarMass                         = array('f', [-99.]);  newbr.append(bBarMass)
bBarMassE                        = array('f', [-99.]);  newbr.append(bBarMassE)
bPx                              = array('f', [-99.]);  newbr.append(bPx)
bPy                              = array('f', [-99.]);  newbr.append(bPy)
bPz                              = array('f', [-99.]);  newbr.append(bPz)
bVtxCL                           = array('f', [-99.]);  newbr.append(bVtxCL)
bVtxX                            = array('f', [-99.]);  newbr.append(bVtxX)
bVtxY                            = array('f', [-99.]);  newbr.append(bVtxY)
bVtxZ                            = array('f', [-99.]);  newbr.append(bVtxZ)
bCosAlphaBS                      = array('f', [-99.]);  newbr.append(bCosAlphaBS)
bCosAlphaBSE                     = array('f', [-99.]);  newbr.append(bCosAlphaBSE)
bLBS                             = array('f', [-99.]);  newbr.append(bLBS)
bLBSE                            = array('f', [-99.]);  newbr.append(bLBSE)
bDCABS                           = array('f', [-99.]);  newbr.append(bDCABS)
bDCABSE                          = array('f', [-99.]);  newbr.append(bDCABSE)
# bctauPVBS                        = array('f', [-99.]);  newbr.append(bctauPVBS)
# bctauPVBSE                       = array('f', [-99.]);  newbr.append(bctauPVBSE)
kstMass                          = array('f', [-99.]);  newbr.append(kstMass)
kstMassE                         = array('f', [-99.]);  newbr.append(kstMassE)
kstBarMass                       = array('f', [-99.]);  newbr.append(kstBarMass)
kstBarMassE                      = array('f', [-99.]);  newbr.append(kstBarMassE)
kstPx                            = array('f', [-99.]);  newbr.append(kstPx)
kstPy                            = array('f', [-99.]);  newbr.append(kstPy)
kstPz                            = array('f', [-99.]);  newbr.append(kstPz)
# kstPxxE                          = array('f', [-99.]);  newbr.append(kstPxxE)
# kstPyyE                          = array('f', [-99.]);  newbr.append(kstPyyE)
# kstPzzE                          = array('f', [-99.]);  newbr.append(kstPzzE)
# kstPxyE                          = array('f', [-99.]);  newbr.append(kstPxyE)
# kstPxzE                          = array('f', [-99.]);  newbr.append(kstPxzE)
# kstPyzE                          = array('f', [-99.]);  newbr.append(kstPyzE)
kstVtxCL                         = array('f', [-99.]);  newbr.append(kstVtxCL)
kstVtxX                          = array('f', [-99.]);  newbr.append(kstVtxX)
kstVtxY                          = array('f', [-99.]);  newbr.append(kstVtxY)
kstVtxZ                          = array('f', [-99.]);  newbr.append(kstVtxZ)
kkMass                           = array('f', [-99.]);  newbr.append(kkMass)
mumuMass                         = array('f', [-99.]);  newbr.append(mumuMass)
mumuMassE                        = array('f', [-99.]);  newbr.append(mumuMassE)
mumuPx                           = array('f', [-99.]);  newbr.append(mumuPx)
mumuPy                           = array('f', [-99.]);  newbr.append(mumuPy)
mumuPz                           = array('f', [-99.]);  newbr.append(mumuPz)
mumuVtxCL                        = array('f', [-99.]);  newbr.append(mumuVtxCL)
mumuVtxX                         = array('f', [-99.]);  newbr.append(mumuVtxX)
mumuVtxY                         = array('f', [-99.]);  newbr.append(mumuVtxY)
mumuVtxZ                         = array('f', [-99.]);  newbr.append(mumuVtxZ)
mumuCosAlphaBS                   = array('f', [-99.]);  newbr.append(mumuCosAlphaBS)
mumuCosAlphaBSE                  = array('f', [-99.]);  newbr.append(mumuCosAlphaBSE)
mumuLBS                          = array('f', [-99.]);  newbr.append(mumuLBS)
mumuLBSE                         = array('f', [-99.]);  newbr.append(mumuLBSE)
mumuDCA                          = array('f', [-99.]);  newbr.append(mumuDCA)
mumHighPurity                    = array('f', [-99.]);  newbr.append(mumHighPurity)
mumCL                            = array('f', [-99.]);  newbr.append(mumCL)
mumNormChi2                      = array('f', [-99.]);  newbr.append(mumNormChi2)
mumPx                            = array('f', [-99.]);  newbr.append(mumPx)
mumPy                            = array('f', [-99.]);  newbr.append(mumPy)
mumPz                            = array('f', [-99.]);  newbr.append(mumPz)
# mumDCAVtx                        = array('f', [-99.]);  newbr.append(mumDCAVtx)
# mumDCAVtxE                       = array('f', [-99.]);  newbr.append(mumDCAVtxE)
mumDCABS                         = array('f', [-99.]);  newbr.append(mumDCABS)
mumDCABSE                        = array('f', [-99.]);  newbr.append(mumDCABSE)
# mumdxyVtx                        = array('f', [-99.]);  newbr.append(mumdxyVtx)
# mumdzVtx                         = array('f', [-99.]);  newbr.append(mumdzVtx)
# mumMinIP2D                       = array('f', [-99.]);  newbr.append(mumMinIP2D)
# mumMinIP2DE                      = array('f', [-99.]);  newbr.append(mumMinIP2DE)
mumNPixLayers                    = array('f', [-99.]);  newbr.append(mumNPixLayers)
mumNTrkLayers                    = array('f', [-99.]);  newbr.append(mumNTrkLayers)
# mumNMuonHits                     = array('f', [-99.]);  newbr.append(mumNMuonHits)
mumNMatchStation                 = array('f', [-99.]);  newbr.append(mumNMatchStation)
mupHighPurity                    = array('f', [-99.]);  newbr.append(mupHighPurity)
mupCL                            = array('f', [-99.]);  newbr.append(mupCL)
mupNormChi2                      = array('f', [-99.]);  newbr.append(mupNormChi2)
mupPx                            = array('f', [-99.]);  newbr.append(mupPx)
mupPy                            = array('f', [-99.]);  newbr.append(mupPy)
mupPz                            = array('f', [-99.]);  newbr.append(mupPz)
# mupDCAVtx                        = array('f', [-99.]);  newbr.append(mupDCAVtx)
# mupDCAVtxE                       = array('f', [-99.]);  newbr.append(mupDCAVtxE)
mupDCABS                         = array('f', [-99.]);  newbr.append(mupDCABS)
mupDCABSE                        = array('f', [-99.]);  newbr.append(mupDCABSE)
mupdxyVtx                        = array('f', [-99.]);  newbr.append(mupdxyVtx)
mupdzVtx                         = array('f', [-99.]);  newbr.append(mupdzVtx)
# mupMinIP2D                       = array('f', [-99.]);  newbr.append(mupMinIP2D)
# mupMinIP2DE                      = array('f', [-99.]);  newbr.append(mupMinIP2DE)
mupNPixLayers                    = array('f', [-99.]);  newbr.append(mupNPixLayers)
mupNTrkLayers                    = array('f', [-99.]);  newbr.append(mupNTrkLayers)
# mupNMuonHits                     = array('f', [-99.]);  newbr.append(mupNMuonHits)
mupNMatchStation                 = array('f', [-99.]);  newbr.append(mupNMatchStation)
kstTrkmHighPurity                = array('f', [-99.]);  newbr.append(kstTrkmHighPurity)
kstTrkmCL                        = array('f', [-99.]);  newbr.append(kstTrkmCL)
kstTrkmNormChi2                  = array('f', [-99.]);  newbr.append(kstTrkmNormChi2)
kstTrkmPx                        = array('f', [-99.]);  newbr.append(kstTrkmPx)
kstTrkmPy                        = array('f', [-99.]);  newbr.append(kstTrkmPy)
kstTrkmPz                        = array('f', [-99.]);  newbr.append(kstTrkmPz)
# kstTrkmPxxE                      = array('f', [-99.]);  newbr.append(kstTrkmPxxE)
# kstTrkmPyyE                      = array('f', [-99.]);  newbr.append(kstTrkmPyyE)
# kstTrkmPzzE                      = array('f', [-99.]);  newbr.append(kstTrkmPzzE)
# kstTrkmPxyE                      = array('f', [-99.]);  newbr.append(kstTrkmPxyE)
# kstTrkmPxzE                      = array('f', [-99.]);  newbr.append(kstTrkmPxzE)
# kstTrkmPyzE                      = array('f', [-99.]);  newbr.append(kstTrkmPyzE)
kstTrkmDCABS                     = array('f', [-99.]);  newbr.append(kstTrkmDCABS)
kstTrkmDCABSE                    = array('f', [-99.]);  newbr.append(kstTrkmDCABSE)
kstTrkmFracHits                  = array('f', [-99.]);  newbr.append(kstTrkmFracHits)
kstTrkmdxyVtx                    = array('f', [-99.]);  newbr.append(kstTrkmdxyVtx)
kstTrkmdzVtx                     = array('f', [-99.]);  newbr.append(kstTrkmdzVtx)
# kstTrkmMinIP2D                   = array('f', [-99.]);  newbr.append(kstTrkmMinIP2D)
# kstTrkmMinIP2DE                  = array('f', [-99.]);  newbr.append(kstTrkmMinIP2DE)
kstTrkmNPixHits                  = array('f', [-99.]);  newbr.append(kstTrkmNPixHits)
kstTrkmNPixLayers                = array('f', [-99.]);  newbr.append(kstTrkmNPixLayers)
kstTrkmNTrkHits                  = array('f', [-99.]);  newbr.append(kstTrkmNTrkHits)
kstTrkmNTrkLayers                = array('f', [-99.]);  newbr.append(kstTrkmNTrkLayers)
kstTrkmMuMatch                   = array('f', [-99.]);  newbr.append(kstTrkmMuMatch)
kstTrkpHighPurity                = array('f', [-99.]);  newbr.append(kstTrkpHighPurity)
kstTrkpCL                        = array('f', [-99.]);  newbr.append(kstTrkpCL)
kstTrkpNormChi2                  = array('f', [-99.]);  newbr.append(kstTrkpNormChi2)
kstTrkpPx                        = array('f', [-99.]);  newbr.append(kstTrkpPx)
kstTrkpPy                        = array('f', [-99.]);  newbr.append(kstTrkpPy)
kstTrkpPz                        = array('f', [-99.]);  newbr.append(kstTrkpPz)
# kstTrkpPxxE                      = array('f', [-99.]);  newbr.append(kstTrkpPxxE)
# kstTrkpPyyE                      = array('f', [-99.]);  newbr.append(kstTrkpPyyE)
# kstTrkpPzzE                      = array('f', [-99.]);  newbr.append(kstTrkpPzzE)
# kstTrkpPxyE                      = array('f', [-99.]);  newbr.append(kstTrkpPxyE)
# kstTrkpPxzE                      = array('f', [-99.]);  newbr.append(kstTrkpPxzE)
# kstTrkpPyzE                      = array('f', [-99.]);  newbr.append(kstTrkpPyzE)
kstTrkpDCABS                     = array('f', [-99.]);  newbr.append(kstTrkpDCABS)
kstTrkpDCABSE                    = array('f', [-99.]);  newbr.append(kstTrkpDCABSE)
kstTrkpFracHits                  = array('f', [-99.]);  newbr.append(kstTrkpFracHits)
kstTrkpdxyVtx                    = array('f', [-99.]);  newbr.append(kstTrkpdxyVtx)
kstTrkpdzVtx                     = array('f', [-99.]);  newbr.append(kstTrkpdzVtx)
# kstTrkpMinIP2D                   = array('f', [-99.]);  newbr.append(kstTrkpMinIP2D)
# kstTrkpMinIP2DE                  = array('f', [-99.]);  newbr.append(kstTrkpMinIP2DE)
kstTrkpNPixHits                  = array('f', [-99.]);  newbr.append(kstTrkpNPixHits)
kstTrkpNPixLayers                = array('f', [-99.]);  newbr.append(kstTrkpNPixLayers)
kstTrkpNTrkHits                  = array('f', [-99.]);  newbr.append(kstTrkpNTrkHits)
kstTrkpNTrkLayers                = array('f', [-99.]);  newbr.append(kstTrkpNTrkLayers)
kstTrkpMuMatch                   = array('f', [-99.]);  newbr.append(kstTrkpMuMatch)
tagB0                            = array('f', [-99.]);  newbr.append(tagB0)

truthMatchMum                    = array('f', [-99.]);  newbr.append(truthMatchMum)
truthMatchMup                    = array('f', [-99.]);  newbr.append(truthMatchMup)
truthMatchTrkm                   = array('f', [-99.]);  newbr.append(truthMatchTrkm)
truthMatchTrkp                   = array('f', [-99.]);  newbr.append(truthMatchTrkp)
mumDeltaRwithMC                  = array('f', [-99.]);  newbr.append(mumDeltaRwithMC)
mupDeltaRwithMC                  = array('f', [-99.]);  newbr.append(mupDeltaRwithMC)
kstTrkpDeltaRwithMC              = array('f', [-99.]);  newbr.append(kstTrkpDeltaRwithMC)
kstTrkmDeltaRwithMC              = array('f', [-99.]);  newbr.append(kstTrkmDeltaRwithMC)
genSignal                        = array('f', [-99.]);  newbr.append(genSignal)
genSignHasFSR                    = array('f', [-99.]);  newbr.append(genSignHasFSR)
genSignKstHasFSR                 = array('f', [-99.]);  newbr.append(genSignKstHasFSR)
genSignPsiHasFSR                 = array('f', [-99.]);  newbr.append(genSignPsiHasFSR)
genPriVtxX                       = array('f', [-99.]);  newbr.append(genPriVtxX)
genPriVtxY                       = array('f', [-99.]);  newbr.append(genPriVtxY)
genPriVtxZ                       = array('f', [-99.]);  newbr.append(genPriVtxZ)
genB0Mass                        = array('f', [-99.]);  newbr.append(genB0Mass)
genB0Px                          = array('f', [-99.]);  newbr.append(genB0Px)
genB0Py                          = array('f', [-99.]);  newbr.append(genB0Py)
genB0Pz                          = array('f', [-99.]);  newbr.append(genB0Pz)
genB0VtxX                        = array('f', [-99.]);  newbr.append(genB0VtxX)
genB0VtxY                        = array('f', [-99.]);  newbr.append(genB0VtxY)
genB0VtxZ                        = array('f', [-99.]);  newbr.append(genB0VtxZ)
genKstMass                       = array('f', [-99.]);  newbr.append(genKstMass)
genKstPx                         = array('f', [-99.]);  newbr.append(genKstPx)
genKstPy                         = array('f', [-99.]);  newbr.append(genKstPy)
genKstPz                         = array('f', [-99.]);  newbr.append(genKstPz)
genKstVtxX                       = array('f', [-99.]);  newbr.append(genKstVtxX)
genKstVtxY                       = array('f', [-99.]);  newbr.append(genKstVtxY)
genKstVtxZ                       = array('f', [-99.]);  newbr.append(genKstVtxZ)
genPsiMass                       = array('f', [-99.]);  newbr.append(genPsiMass)
genPsiVtxX                       = array('f', [-99.]);  newbr.append(genPsiVtxX)
genPsiVtxY                       = array('f', [-99.]);  newbr.append(genPsiVtxY)
genPsiVtxZ                       = array('f', [-99.]);  newbr.append(genPsiVtxZ)
genMumMother                     = array('f', [-99.]);  newbr.append(genMumMother)
genMumPx                         = array('f', [-99.]);  newbr.append(genMumPx)
genMumPy                         = array('f', [-99.]);  newbr.append(genMumPy)
genMumPz                         = array('f', [-99.]);  newbr.append(genMumPz)
genMupMother                     = array('f', [-99.]);  newbr.append(genMupMother)
genMupPx                         = array('f', [-99.]);  newbr.append(genMupPx)
genMupPy                         = array('f', [-99.]);  newbr.append(genMupPy)
genMupPz                         = array('f', [-99.]);  newbr.append(genMupPz)
genKstTrkmMother                 = array('f', [-99.]);  newbr.append(genKstTrkmMother)
genKstTrkmID                     = array('f', [-99.]);  newbr.append(genKstTrkmID)
genKstTrkmPx                     = array('f', [-99.]);  newbr.append(genKstTrkmPx)
genKstTrkmPy                     = array('f', [-99.]);  newbr.append(genKstTrkmPy)
genKstTrkmPz                     = array('f', [-99.]);  newbr.append(genKstTrkmPz)
genKstTrkpMother                 = array('f', [-99.]);  newbr.append(genKstTrkpMother)
genKstTrkpID                     = array('f', [-99.]);  newbr.append(genKstTrkpID)
genKstTrkpPx                     = array('f', [-99.]);  newbr.append(genKstTrkpPx)
genKstTrkpPy                     = array('f', [-99.]);  newbr.append(genKstTrkpPy)
genKstTrkpPz                     = array('f', [-99.]);  newbr.append(genKstTrkpPz)

mumGlobalMuon                    = array('f',[-99.]); newbr.append(mumGlobalMuon             )   
mumTrackerMuon                   = array('f',[-99.]); newbr.append(mumTrackerMuon            )   
mumStandAloneMuon                = array('f',[-99.]); newbr.append(mumStandAloneMuon         )   
mumTMOneStationTight             = array('f',[-99.]); newbr.append(mumTMOneStationTight      )   
mumTMOneStationLoose             = array('f',[-99.]); newbr.append(mumTMOneStationLoose      )   
mupGlobalMuon                    = array('f',[-99.]); newbr.append(mupGlobalMuon             )   
mupTrackerMuon                   = array('f',[-99.]); newbr.append(mupTrackerMuon            )   
mupStandAloneMuon                = array('f',[-99.]); newbr.append(mupStandAloneMuon         )   
mupTMOneStationTight             = array('f',[-99.]); newbr.append(mupTMOneStationTight      )   
mupTMOneStationLoose             = array('f',[-99.]); newbr.append(mupTMOneStationLoose      )   
kstTrkmGlobalMuon                = array('f',[-99.]); newbr.append(kstTrkmGlobalMuon             )   
kstTrkmTrackerMuon               = array('f',[-99.]); newbr.append(kstTrkmTrackerMuon            )   
kstTrkmStandAloneMuon            = array('f',[-99.]); newbr.append(kstTrkmStandAloneMuon         )   
kstTrkmGlobalMuonPromptTight     = array('f',[-99.]); newbr.append(kstTrkmGlobalMuonPromptTight  )   
kstTrkmTMOneStationTight         = array('f',[-99.]); newbr.append(kstTrkmTMOneStationTight      )   
kstTrkmTMOneStationLoose         = array('f',[-99.]); newbr.append(kstTrkmTMOneStationLoose      )   
kstTrkmTrackerMuonArbitrated     = array('f',[-99.]); newbr.append(kstTrkmTrackerMuonArbitrated  )   
kstTrkpGlobalMuon                = array('f',[-99.]); newbr.append(kstTrkpGlobalMuon             )   
kstTrkpTrackerMuon               = array('f',[-99.]); newbr.append(kstTrkpTrackerMuon            )   
kstTrkpStandAloneMuon            = array('f',[-99.]); newbr.append(kstTrkpStandAloneMuon         )   
kstTrkpGlobalMuonPromptTight     = array('f',[-99.]); newbr.append(kstTrkpGlobalMuonPromptTight  )   
kstTrkpTMOneStationTight         = array('f',[-99.]); newbr.append(kstTrkpTMOneStationTight      )   
kstTrkpTMOneStationLoose         = array('f',[-99.]); newbr.append(kstTrkpTMOneStationLoose      )   
kstTrkpTrackerMuonArbitrated     = array('f',[-99.]); newbr.append(kstTrkpTrackerMuonArbitrated  )   


bPt                              = array('f',[-99.]); newbr.append( bPt          )   
kstPt                            = array('f',[-99.]); newbr.append( kstPt        )   
mumuPt                           = array('f',[-99.]); newbr.append( mumuPt       )   
mumPt                            = array('f',[-99.]); newbr.append( mumPt        )   
mupPt                            = array('f',[-99.]); newbr.append( mupPt        )   
kstTrkmPt                        = array('f',[-99.]); newbr.append( kstTrkmPt    )   
kstTrkpPt                        = array('f',[-99.]); newbr.append( kstTrkpPt    )   
bPhi                             = array('f',[-99.]); newbr.append( bPhi         )   
kstPhi                           = array('f',[-99.]); newbr.append( kstPhi       )   
mumuPhi                          = array('f',[-99.]); newbr.append( mumuPhi      )   
mumPhi                           = array('f',[-99.]); newbr.append( mumPhi       )   
mupPhi                           = array('f',[-99.]); newbr.append( mupPhi       )   
kstTrkmPhi                       = array('f',[-99.]); newbr.append( kstTrkmPhi   )   
kstTrkpPhi                       = array('f',[-99.]); newbr.append( kstTrkpPhi   )   
bEta                             = array('f',[-99.]); newbr.append( bEta         )   
kstEta                           = array('f',[-99.]); newbr.append( kstEta       )   
mumuEta                          = array('f',[-99.]); newbr.append( mumuEta      )   
mumEta                           = array('f',[-99.]); newbr.append( mumEta       )   
mupEta                           = array('f',[-99.]); newbr.append( mupEta       )   
kstTrkmEta                       = array('f',[-99.]); newbr.append( kstTrkmEta   )   
kstTrkpEta                       = array('f',[-99.]); newbr.append( kstTrkpEta   )   


# mumIsoN                          = ROOT.std.vector('float')()   
# mupIsoN                          = ROOT.std.vector('float')()   
# kstTrkmIsoN                      = ROOT.std.vector('float')()   
# kstTrkpIsoN                      = ROOT.std.vector('float')()   

if not skim:
  bCosAlphaVtx                     = array('f', [-99.]);  newbr.append(bCosAlphaVtx)
  bCosAlphaVtxE                    = array('f', [-99.]);  newbr.append(bCosAlphaVtxE)
  priVtxCL                         = array('f', [-99.]);  newbr.append(priVtxCL)
  priVtxX                          = array('f', [-99.]);  newbr.append(priVtxX)
  priVtxY                          = array('f', [-99.]);  newbr.append(priVtxY)
  priVtxZ                          = array('f', [-99.]);  newbr.append(priVtxZ)
#   bLVtx                            = array('f', [-99.]);  newbr.append(bLVtx)
#   bLVtxE                           = array('f', [-99.]);  newbr.append(bLVtxE)
#   bDCAVtx                          = array('f', [-99.]);  newbr.append(bDCAVtx)
#   bDCAVtxE                         = array('f', [-99.]);  newbr.append(bDCAVtxE)
  mumNPixHits                      = array('f', [-99.]);  newbr.append(mumNPixHits)
  mumNTrkHits                      = array('f', [-99.]);  newbr.append(mumNTrkHits)
  mupNPixHits                      = array('f', [-99.]);  newbr.append(mupNPixHits)
  mupNTrkHits                      = array('f', [-99.]);  newbr.append(mupNTrkHits)
#   mumKinkChi2                      = array('f', [-99.]);  newbr.append(mumKinkChi2)
  mumFracHits                      = array('f', [-99.]);  newbr.append(mumFracHits)
#   mupKinkChi2                      = array('f', [-99.]);  newbr.append(mupKinkChi2)
  mupFracHits                      = array('f', [-99.]);  newbr.append(mupFracHits)
#   mumMinIP                         = array('f', [-99.]);  newbr.append(mumMinIP)
#   mumMinIPE                        = array('f', [-99.]);  newbr.append(mumMinIPE)
#   mupMinIP                         = array('f', [-99.]);  newbr.append(mupMinIP)
#   mupMinIPE                        = array('f', [-99.]);  newbr.append(mupMinIPE)
#   kstTrkmDCAVtx                    = array('f', [-99.]);  newbr.append(kstTrkmDCAVtx)
#   kstTrkmDCAVtxE                   = array('f', [-99.]);  newbr.append(kstTrkmDCAVtxE)
#   kstTrkpDCAVtx                    = array('f', [-99.]);  newbr.append(kstTrkpDCAVtx)
#   kstTrkpDCAVtxE                   = array('f', [-99.]);  newbr.append(kstTrkpDCAVtxE)
  mumGlobalMuonPromptTight         = array('f',[-99.]); newbr.append(mumGlobalMuonPromptTight  )   
  mumTMLastStationTight            = array('f',[-99.]); newbr.append(mumTMLastStationTight     )   
  mumTMLastStationLoose            = array('f',[-99.]); newbr.append(mumTMLastStationLoose     )   
  mumTM2DCompatibilityTight        = array('f',[-99.]); newbr.append(mumTM2DCompatibilityTight )   
  mumTM2DCompatibilityLoose        = array('f',[-99.]); newbr.append(mumTM2DCompatibilityLoose )   
  mumTMLastStationAngTight         = array('f',[-99.]); newbr.append(mumTMLastStationAngTight  )   
  mumTMLastStationAngLoose         = array('f',[-99.]); newbr.append(mumTMLastStationAngLoose  )   
  mumTMOneStationAngTight          = array('f',[-99.]); newbr.append(mumTMOneStationAngTight   )   
  mumTMOneStationAngLoose          = array('f',[-99.]); newbr.append(mumTMOneStationAngLoose   )   
  mumTrackerMuonArbitrated         = array('f',[-99.]); newbr.append(mumTrackerMuonArbitrated  )   
  mupGlobalMuonPromptTight         = array('f',[-99.]); newbr.append(mupGlobalMuonPromptTight  )   
  mupTMLastStationTight            = array('f',[-99.]); newbr.append(mupTMLastStationTight     )   
  mupTMLastStationLoose            = array('f',[-99.]); newbr.append(mupTMLastStationLoose     )   
  mupTM2DCompatibilityTight        = array('f',[-99.]); newbr.append(mupTM2DCompatibilityTight )   
  mupTM2DCompatibilityLoose        = array('f',[-99.]); newbr.append(mupTM2DCompatibilityLoose )   
  mupTMLastStationAngTight         = array('f',[-99.]); newbr.append(mupTMLastStationAngTight  )   
  mupTMLastStationAngLoose         = array('f',[-99.]); newbr.append(mupTMLastStationAngLoose  )   
  mupTMOneStationAngTight          = array('f',[-99.]); newbr.append(mupTMOneStationAngTight   )   
  mupTMOneStationAngLoose          = array('f',[-99.]); newbr.append(mupTMOneStationAngLoose   )   
  mupTrackerMuonArbitrated         = array('f',[-99.]); newbr.append(mupTrackerMuonArbitrated  )   
  kstTrkmTM2DCompatibilityTight    = array('f',[-99.]); newbr.append(kstTrkmTM2DCompatibilityTight )   
  kstTrkmTM2DCompatibilityLoose    = array('f',[-99.]); newbr.append(kstTrkmTM2DCompatibilityLoose )   
  kstTrkpTM2DCompatibilityTight    = array('f',[-99.]); newbr.append(kstTrkpTM2DCompatibilityTight )   
  kstTrkpTM2DCompatibilityLoose    = array('f',[-99.]); newbr.append(kstTrkpTM2DCompatibilityLoose )   
  kstTrkmTMLastStationAngTight     = array('f',[-99.]); newbr.append(kstTrkmTMLastStationAngTight  )   
  kstTrkmTMLastStationAngLoose     = array('f',[-99.]); newbr.append(kstTrkmTMLastStationAngLoose  )   
  kstTrkmTMOneStationAngTight      = array('f',[-99.]); newbr.append(kstTrkmTMOneStationAngTight   )   
  kstTrkmTMOneStationAngLoose      = array('f',[-99.]); newbr.append(kstTrkmTMOneStationAngLoose   )   
  kstTrkpTMLastStationAngTight     = array('f',[-99.]); newbr.append(kstTrkpTMLastStationAngTight  )   
  kstTrkpTMLastStationAngLoose     = array('f',[-99.]); newbr.append(kstTrkpTMLastStationAngLoose  )   
  kstTrkpTMOneStationAngTight      = array('f',[-99.]); newbr.append(kstTrkpTMOneStationAngTight   )   
  kstTrkpTMOneStationAngLoose      = array('f',[-99.]); newbr.append(kstTrkpTMOneStationAngLoose   )   
  kstTrkmTMLastStationTight        = array('f',[-99.]); newbr.append(kstTrkmTMLastStationTight     )   
  kstTrkmTMLastStationLoose        = array('f',[-99.]); newbr.append(kstTrkmTMLastStationLoose     )   
  kstTrkpTMLastStationTight        = array('f',[-99.]); newbr.append(kstTrkpTMLastStationTight     )   
  kstTrkpTMLastStationLoose        = array('f',[-99.]); newbr.append(kstTrkpTMLastStationLoose     )   
#   kstTrkmMinIP                     = array('f', [-99.]);  newbr.append(kstTrkmMinIP)
#   kstTrkmMinIPE                    = array('f', [-99.]);  newbr.append(kstTrkmMinIPE)
#   kstTrkpMinIP                     = array('f', [-99.]);  newbr.append(kstTrkpMinIP)
#   kstTrkpMinIPE                    = array('f', [-99.]);  newbr.append(kstTrkpMinIPE)

#   mumIso                           = ROOT.std.vector('float')()  
#   mupIso                           = ROOT.std.vector('float')()  
#   kstTrkmIso                       = ROOT.std.vector('float')()  
#   kstTrkpIso                       = ROOT.std.vector('float')()  


'''
for i in branches:
    print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/F')".format(NAME=i).expandtabs(40)
'''
ntuple.Branch('runN',                   runN,                                   'runN/F')
ntuple.Branch('eventN',                 eventN,                                 'eventN/F')
ntuple.Branch('recoVtxN',               recoVtxN,                               'recoVtxN/F')
ntuple.Branch('evWeight',               evWeight,                               'evWeight/F')
ntuple.Branch('evWeightE2',             evWeightE2,                             'evWeightE2/F')
ntuple.Branch('trueNumInteractionsMC',  trueNumInteractionsMC,                  'trueNumInteractionsMC/F')
ntuple.Branch('trig',                   trig,                                   'trig/F')
ntuple.Branch('bsX',                    bsX,                                    'bsX/F')
ntuple.Branch('bsY',                    bsY,                                    'bsY/F')
ntuple.Branch('bMass',                  bMass,                                  'bMass/F')
ntuple.Branch('bMassE',                 bMassE,                                 'bMassE/F')
ntuple.Branch('bBarMass',               bBarMass,                               'bBarMass/F')
ntuple.Branch('bBarMassE',              bBarMassE,                              'bBarMassE/F')
ntuple.Branch('bPx',                    bPx,                                    'bPx/F')
ntuple.Branch('bPy',                    bPy,                                    'bPy/F')
ntuple.Branch('bPz',                    bPz,                                    'bPz/F')
ntuple.Branch('bVtxCL',                 bVtxCL,                                 'bVtxCL/F')
ntuple.Branch('bVtxX',                  bVtxX,                                  'bVtxX/F')
ntuple.Branch('bVtxY',                  bVtxY,                                  'bVtxY/F')
ntuple.Branch('bVtxZ',                  bVtxZ,                                  'bVtxZ/F')

ntuple.Branch('bCosAlphaBS',            bCosAlphaBS,                            'bCosAlphaBS/F')
ntuple.Branch('bCosAlphaBSE',           bCosAlphaBSE,                           'bCosAlphaBSE/F')
ntuple.Branch('bLBS',                   bLBS,                                   'bLBS/F')
ntuple.Branch('bLBSE',                  bLBSE,                                  'bLBSE/F')
ntuple.Branch('bDCABS',                 bDCABS,                                 'bDCABS/F')
ntuple.Branch('bDCABSE',                bDCABSE,                                'bDCABSE/F')
# ntuple.Branch('bctauPVBS',              bctauPVBS,                              'bctauPVBS/F')
# ntuple.Branch('bctauPVBSE',             bctauPVBSE,                             'bctauPVBSE/F')
ntuple.Branch('kstMass',                kstMass,                                'kstMass/F')
ntuple.Branch('kstMassE',               kstMassE,                               'kstMassE/F')
ntuple.Branch('kstBarMass',             kstBarMass,                             'kstBarMass/F')
ntuple.Branch('kstBarMassE',            kstBarMassE,                            'kstBarMassE/F')
ntuple.Branch('kstPx',                  kstPx,                                  'kstPx/F')
ntuple.Branch('kstPy',                  kstPy,                                  'kstPy/F')
ntuple.Branch('kstPz',                  kstPz,                                  'kstPz/F')
# ntuple.Branch('kstPxxE',                kstPxxE,                                'kstPxxE/F')
# ntuple.Branch('kstPyyE',                kstPyyE,                                'kstPyyE/F')
# ntuple.Branch('kstPzzE',                kstPzzE,                                'kstPzzE/F')
# ntuple.Branch('kstPxyE',                kstPxyE,                                'kstPxyE/F')
# ntuple.Branch('kstPxzE',                kstPxzE,                                'kstPxzE/F')
# ntuple.Branch('kstPyzE',                kstPyzE,                                'kstPyzE/F')
ntuple.Branch('kstVtxCL',               kstVtxCL,                               'kstVCL/F')
ntuple.Branch('kstVtxX',                kstVtxX,                                'kstVtxX/F')
ntuple.Branch('kstVtxY',                kstVtxY,                                'kstVtxY/F')
ntuple.Branch('kstVtxZ',                kstVtxZ,                                'kstVtxZ/F')
ntuple.Branch('kkMass',                 kkMass,                                 'kkMass/F')
ntuple.Branch('mumuMass',               mumuMass,                               'mumuMass/F')
ntuple.Branch('mumuMassE',              mumuMassE,                              'mumuMassE/F')
ntuple.Branch('mumuPx',                 mumuPx,                                 'mumuPx/F')
ntuple.Branch('mumuPy',                 mumuPy,                                 'mumuPy/F')
ntuple.Branch('mumuPz',                 mumuPz,                                 'mumuPz/F')
ntuple.Branch('mumuVtxCL',              mumuVtxCL,                              'mumuVtxCL/F')
ntuple.Branch('mumuVtxX',               mumuVtxX,                               'mumuVtxX/F')
ntuple.Branch('mumuVtxY',               mumuVtxY,                               'mumuVtxY/F')
ntuple.Branch('mumuVtxZ',               mumuVtxZ,                               'mumuVtxZ/F')
ntuple.Branch('mumuCosAlphaBS',         mumuCosAlphaBS,                         'mumuCosAlphaBS/F')
ntuple.Branch('mumuCosAlphaBSE',        mumuCosAlphaBSE,                        'mumuCosAlphaBSE/F')
ntuple.Branch('mumuLBS',                mumuLBS,                                'mumuLBS/F')
ntuple.Branch('mumuLBSE',               mumuLBSE,                               'mumuLBSE/F')
ntuple.Branch('mumuDCA',                mumuDCA,                                'mumuDCA/F')

ntuple.Branch('mumHighPurity',          mumHighPurity,                          'mumHighPurity/F')
ntuple.Branch('mumCL',                  mumCL,                                  'mumCL/F')
ntuple.Branch('mumNormChi2',            mumNormChi2,                            'mumNormChi2/F')
ntuple.Branch('mumPx',                  mumPx,                                  'mumPx/F')
ntuple.Branch('mumPy',                  mumPy,                                  'mumPy/F')
ntuple.Branch('mumPz',                  mumPz,                                  'mumPz/F')
# ntuple.Branch('mumDCAVtx',              mumDCAVtx,                              'mumDCAVtx/F')
# ntuple.Branch('mumDCAVtxE',             mumDCAVtxE,                             'mumDCAVtxE/F')
ntuple.Branch('mumDCABS',               mumDCABS,                               'mumDCABS/F')
ntuple.Branch('mumDCABSE',              mumDCABSE,                              'mumDCABSE/F')
# ntuple.Branch('mumdxyVtx',              mumdxyVtx,                              'mumdxyVtx/F')
# ntuple.Branch('mumdzVtx',               mumdzVtx,                               'mumdzVtx/F')
# ntuple.Branch('mumMinIP2D',             mumMinIP2D,                             'mumMinIP2D/F')
# ntuple.Branch('mumMinIP2DE',            mumMinIP2DE,                            'mumMinIP2DE/F')
ntuple.Branch('mumNPixLayers',          mumNPixLayers,                          'mumNPixLayers/F')
ntuple.Branch('mumNTrkLayers',          mumNTrkLayers,                          'mumNTrkLayers/F')
# ntuple.Branch('mumNMuonHits',           mumNMuonHits,                           'mumNMuonHits/F')
ntuple.Branch('mumNMatchStation',       mumNMatchStation,                       'mumNMatchStation/F')
ntuple.Branch('mupHighPurity',          mupHighPurity,                          'mupHighPurity/F')
ntuple.Branch('mupCL',                  mupCL,                                  'mupCL/F')
ntuple.Branch('mupNormChi2',            mupNormChi2,                            'mupNormChi2/F')
ntuple.Branch('mupPx',                  mupPx,                                  'mupPx/F')
ntuple.Branch('mupPy',                  mupPy,                                  'mupPy/F')
ntuple.Branch('mupPz',                  mupPz,                                  'mupPz/F')
# ntuple.Branch('mupDCAVtx',              mupDCAVtx,                              'mupDCAVtx/F')
# ntuple.Branch('mupDCAVtxE',             mupDCAVtxE,                             'mupDCAVtxE/F')
ntuple.Branch('mupDCABS',               mupDCABS,                               'mupDCABS/F')
ntuple.Branch('mupDCABSE',              mupDCABSE,                              'mupDCABSE/F')
ntuple.Branch('mupdxyVtx',              mupdxyVtx,                              'mupdxyVtx/F')
ntuple.Branch('mupdzVtx',               mupdzVtx,                               'mupdzVtx/F')
# ntuple.Branch('mupMinIP2D',             mupMinIP2D,                             'mupMinIP2D/F')
# ntuple.Branch('mupMinIP2DE',            mupMinIP2DE,                            'mupMinIP2DE/F')
ntuple.Branch('mupNPixLayers',          mupNPixLayers,                          'mupNPixLayers/F')
ntuple.Branch('mupNTrkLayers',          mupNTrkLayers,                          'mupNTrkLayers/F')
# ntuple.Branch('mupNMuonHits',           mupNMuonHits,                           'mupNMuonHits/F')
ntuple.Branch('mupNMatchStation',       mupNMatchStation,                       'mupNMatchStation/F')

ntuple.Branch('kstTrkmHighPurity',      kstTrkmHighPurity,                      'kstTrkmHighPurity/F')
ntuple.Branch('kstTrkmCL',              kstTrkmCL,                              'kstTrkmCL/F')
ntuple.Branch('kstTrkmNormChi2',        kstTrkmNormChi2,                        'kstTrkmNormChi2/F')
ntuple.Branch('kstTrkmPx',              kstTrkmPx,                              'kstTrkmPx/F')
ntuple.Branch('kstTrkmPy',              kstTrkmPy,                              'kstTrkmPy/F')
ntuple.Branch('kstTrkmPz',              kstTrkmPz,                              'kstTrkmPz/F')
# ntuple.Branch('kstTrkmPxxE',            kstTrkmPxxE,                            'kstTrkmPxxE/F')
# ntuple.Branch('kstTrkmPyyE',            kstTrkmPyyE,                            'kstTrkmPyyE/F')
# ntuple.Branch('kstTrkmPzzE',            kstTrkmPzzE,                            'kstTrkmPzzE/F')
# ntuple.Branch('kstTrkmPxyE',            kstTrkmPxyE,                            'kstTrkmPxyE/F')
# ntuple.Branch('kstTrkmPxzE',            kstTrkmPxzE,                            'kstTrkmPxzE/F')
# ntuple.Branch('kstTrkmPyzE',            kstTrkmPyzE,                            'kstTrkmPyzE/F')
ntuple.Branch('kstTrkmDCABS',           kstTrkmDCABS,                           'kstTrkmDCABS/F')
ntuple.Branch('kstTrkmDCABSE',          kstTrkmDCABSE,                          'kstTrkmDCABSE/F')
ntuple.Branch('kstTrkmFracHits',        kstTrkmFracHits,                        'kstTrkmFracHits/F')
ntuple.Branch('kstTrkmdxyVtx',          kstTrkmdxyVtx,                          'kstTrkmdxyVtx/F')
ntuple.Branch('kstTrkmdzVtx',           kstTrkmdzVtx,                           'kstTrkmdzVtx/F')
# ntuple.Branch('kstTrkmMinIP2D',         kstTrkmMinIP2D,                         'kstTrkmMinIP2D/F')
# ntuple.Branch('kstTrkmMinIP2DE',        kstTrkmMinIP2DE,                        'kstTrkmMinIP2DE/F')
ntuple.Branch('kstTrkmNPixHits',        kstTrkmNPixHits,                        'kstTrkmNPixHits/F')
ntuple.Branch('kstTrkmNPixLayers',      kstTrkmNPixLayers,                      'kstTrkmNPixLayers/F')
ntuple.Branch('kstTrkmNTrkHits',        kstTrkmNTrkHits,                        'kstTrkmNTrkHits/F')
ntuple.Branch('kstTrkmNTrkLayers',      kstTrkmNTrkLayers,                      'kstTrkmNTrkLayers/F')
ntuple.Branch('kstTrkmMuMatch',         kstTrkmMuMatch,                         'kstTrkmMuMatch/F')
ntuple.Branch('kstTrkpHighPurity',      kstTrkpHighPurity,                      'kstTrkpHighPurity/F')
ntuple.Branch('kstTrkpCL',              kstTrkpCL,                              'kstTrkpCL/F')
ntuple.Branch('kstTrkpNormChi2',        kstTrkpNormChi2,                        'kstTrkpNormChi2/F')
ntuple.Branch('kstTrkpPx',              kstTrkpPx,                              'kstTrkpPx/F')
ntuple.Branch('kstTrkpPy',              kstTrkpPy,                              'kstTrkpPy/F')
ntuple.Branch('kstTrkpPz',              kstTrkpPz,                              'kstTrkpPz/F')
# ntuple.Branch('kstTrkpPxxE',            kstTrkpPxxE,                            'kstTrkpPxxE/F')
# ntuple.Branch('kstTrkpPyyE',            kstTrkpPyyE,                            'kstTrkpPyyE/F')
# ntuple.Branch('kstTrkpPzzE',            kstTrkpPzzE,                            'kstTrkpPzzE/F')
# ntuple.Branch('kstTrkpPxyE',            kstTrkpPxyE,                            'kstTrkpPxyE/F')
# ntuple.Branch('kstTrkpPxzE',            kstTrkpPxzE,                            'kstTrkpPxzE/F')
# ntuple.Branch('kstTrkpPyzE',            kstTrkpPyzE,                            'kstTrkpPyzE/F')
ntuple.Branch('kstTrkpDCABS',           kstTrkpDCABS,                           'kstTrkpDCABS/F')
ntuple.Branch('kstTrkpDCABSE',          kstTrkpDCABSE,                          'kstTrkpDCABSE/F')
ntuple.Branch('kstTrkpFracHits',        kstTrkpFracHits,                        'kstTrkpFracHits/F')
ntuple.Branch('kstTrkpdxyVtx',          kstTrkpdxyVtx,                          'kstTrkpdxyVtx/F')
ntuple.Branch('kstTrkpdzVtx',           kstTrkpdzVtx,                           'kstTrkpdzVtx/F')
# ntuple.Branch('kstTrkpMinIP2D',         kstTrkpMinIP2D,                         'kstTrkpMinIP2D/F')
# ntuple.Branch('kstTrkpMinIP2DE',        kstTrkpMinIP2DE,                        'kstTrkpMinIP2DE/F')
ntuple.Branch('kstTrkpNPixHits',        kstTrkpNPixHits,                        'kstTrkpNPixHits/F')
ntuple.Branch('kstTrkpNPixLayers',      kstTrkpNPixLayers,                      'kstTrkpNPixLayers/F')
ntuple.Branch('kstTrkpNTrkHits',        kstTrkpNTrkHits,                        'kstTrkpNTrkHits/F')
ntuple.Branch('kstTrkpNTrkLayers',      kstTrkpNTrkLayers,                      'kstTrkpNTrkLayers/F')
ntuple.Branch('kstTrkpMuMatch',         kstTrkpMuMatch,                         'kstTrkpMuMatch/F')
ntuple.Branch('tagB0',                  tagB0,                                  'tagB0/F')
ntuple.Branch('truthMatchMum',          truthMatchMum,                          'truthMatchMum/F')
ntuple.Branch('truthMatchMup',          truthMatchMup,                          'truthMatchMup/F')
ntuple.Branch('truthMatchTrkm',         truthMatchTrkm,                         'truthMatchTrkm/F')
ntuple.Branch('truthMatchTrkp',         truthMatchTrkp,                         'truthMatchTrkp/F')
ntuple.Branch('mumDeltaRwithMC',        mumDeltaRwithMC,                        'mumDeltaRwithMC/F')
ntuple.Branch('mupDeltaRwithMC',        mupDeltaRwithMC,                        'mupDeltaRwithMC/F')
ntuple.Branch('kstTrkpDeltaRwithMC',    kstTrkpDeltaRwithMC,                    'kstTrkpDeltaRwithMC/F')
ntuple.Branch('kstTrkmDeltaRwithMC',    kstTrkmDeltaRwithMC,                    'kstTrkmDeltaRwithMC/F')
ntuple.Branch('genSignal',              genSignal,                              'genSignal/F')
ntuple.Branch('genSignHasFSR',          genSignHasFSR,                          'genSignHasFSR/F')
ntuple.Branch('genSignKstHasFSR',       genSignKstHasFSR,                       'genSignKstHasFSR/F')
ntuple.Branch('genSignPsiHasFSR',       genSignPsiHasFSR,                       'genSignPsiHasFSR/F')
ntuple.Branch('genPriVtxX',             genPriVtxX,                             'genPriVtxX/F')
ntuple.Branch('genPriVtxY',             genPriVtxY,                             'genPriVtxY/F')
ntuple.Branch('genPriVtxZ',             genPriVtxZ,                             'genPriVtxZ/F')
ntuple.Branch('genB0Mass',              genB0Mass,                              'genB0Mass/F')
ntuple.Branch('genB0Px',                genB0Px,                                'genB0Px/F')
ntuple.Branch('genB0Py',                genB0Py,                                'genB0Py/F')
ntuple.Branch('genB0Pz',                genB0Pz,                                'genB0Pz/F')
ntuple.Branch('genB0VtxX',              genB0VtxX,                              'genB0VtxX/F')
ntuple.Branch('genB0VtxY',              genB0VtxY,                              'genB0VtxY/F')
ntuple.Branch('genB0VtxZ',              genB0VtxZ,                              'genB0VtxZ/F')
ntuple.Branch('genKstMass',             genKstMass,                             'genKstMass/F')
ntuple.Branch('genKstPx',               genKstPx,                               'genKstPx/F')
ntuple.Branch('genKstPy',               genKstPy,                               'genKstPy/F')
ntuple.Branch('genKstPz',               genKstPz,                               'genKstPz/F')
ntuple.Branch('genKstVtxY',             genKstVtxY,                             'genKstVtxY/F')
ntuple.Branch('genPsiMass',             genPsiMass,                             'genPsiMass/F')
ntuple.Branch('genPsiVtxX',             genPsiVtxX,                             'genPsiVtxX/F')
ntuple.Branch('genPsiVtxY',             genPsiVtxY,                             'genPsiVtxY/F')
ntuple.Branch('genPsiVtxZ',             genPsiVtxZ,                             'genPsiVtxZ/F')
ntuple.Branch('genMumMother',           genMumMother,                           'genMumMother/F')
ntuple.Branch('genMumPx',               genMumPx,                               'genMumPx/F')
ntuple.Branch('genMumPy',               genMumPy,                               'genMumPy/F')
ntuple.Branch('genMumPz',               genMumPz,                               'genMumPz/F')
ntuple.Branch('genMupMother',           genMupMother,                           'genMupMother/F')
ntuple.Branch('genMupPx',               genMupPx,                               'genMupPx/F')
ntuple.Branch('genMupPy',               genMupPy,                               'genMupPy/F')
ntuple.Branch('genMupPz',               genMupPz,                               'genMupPz/F')
ntuple.Branch('genKstTrkmMother',       genKstTrkmMother,                       'genKstTrkmMother/F')
ntuple.Branch('genKstTrkmPx',           genKstTrkmPx,                           'genKstTrkmPx/F')
ntuple.Branch('genKstTrkmPy',           genKstTrkmPy,                           'genKstTrkmPy/F')
ntuple.Branch('genKstTrkmPz',           genKstTrkmPz,                           'genKstTrkmPz/F')
ntuple.Branch('genKstTrkmID',           genKstTrkmID,                           'genKstTrkmID/F')
ntuple.Branch('genKstTrkpMother',       genKstTrkpMother,                       'genKstTrkpMother/F')
ntuple.Branch('genKstTrkpPx',           genKstTrkpPx,                           'genKstTrkpPx/F')
ntuple.Branch('genKstTrkpPy',           genKstTrkpPy,                           'genKstTrkpPy/F')
ntuple.Branch('genKstTrkpPz',           genKstTrkpPz,                           'genKstTrkpPz/F')
ntuple.Branch('genKstTrkpID',           genKstTrkpID,                           'genKstTrkpID/F')

ntuple.Branch('mumGlobalMuon',             mumGlobalMuon              ,               'mumGlobalMuon/F')
ntuple.Branch('mumTrackerMuon',            mumTrackerMuon             ,               'mumTrackerMuon/F')
ntuple.Branch('mumStandAloneMuon',         mumStandAloneMuon          ,               'mumStandAloneMuon/F')
ntuple.Branch('mumTMOneStationTight',      mumTMOneStationTight       ,               'mumTMOneStationTight/F')
ntuple.Branch('mumTMOneStationLoose',      mumTMOneStationLoose       ,               'mumTMOneStationLoose/F')
ntuple.Branch('mupGlobalMuon',             mupGlobalMuon              ,               'mupGlobalMuon/F')
ntuple.Branch('mupTrackerMuon',            mupTrackerMuon             ,               'mupTrackerMuon/F')
ntuple.Branch('mupStandAloneMuon',         mupStandAloneMuon          ,               'mupStandAloneMuon/F')
ntuple.Branch('mupTMOneStationTight',      mupTMOneStationTight       ,               'mupTMOneStationTight/F')
ntuple.Branch('mupTMOneStationLoose',      mupTMOneStationLoose       ,               'mupTMOneStationLoose/F')

ntuple.Branch('kstTrkmGlobalMuon',            kstTrkmGlobalMuon             ,               'kstTrkmGlobalMuon/F')
ntuple.Branch('kstTrkmTrackerMuon',           kstTrkmTrackerMuon            ,               'kstTrkmTrackerMuon/F')
ntuple.Branch('kstTrkmStandAloneMuon',        kstTrkmStandAloneMuon         ,               'kstTrkmStandAloneMuon/F')
ntuple.Branch('kstTrkmGlobalMuonPromptTight', kstTrkmGlobalMuonPromptTight  ,               'kstTrkmGlobalMuonPromptTight/F')
ntuple.Branch('kstTrkmTMOneStationTight',     kstTrkmTMOneStationTight      ,               'kstTrkmTMOneStationTight/F')
ntuple.Branch('kstTrkmTMOneStationLoose',     kstTrkmTMOneStationLoose      ,               'kstTrkmTMOneStationLoose/F')
ntuple.Branch('kstTrkmTrackerMuonArbitrated', kstTrkmTrackerMuonArbitrated  ,               'kstTrkmTrackerMuonArbitrated/F')
ntuple.Branch('kstTrkpGlobalMuon',            kstTrkpGlobalMuon             ,               'kstTrkpGlobalMuon/F')
ntuple.Branch('kstTrkpTrackerMuon',           kstTrkpTrackerMuon            ,               'kstTrkpTrackerMuon/F')
ntuple.Branch('kstTrkpStandAloneMuon',        kstTrkpStandAloneMuon         ,               'kstTrkpStandAloneMuon/F')
ntuple.Branch('kstTrkpGlobalMuonPromptTight', kstTrkpGlobalMuonPromptTight  ,               'kstTrkpGlobalMuonPromptTight/F')
ntuple.Branch('kstTrkpTMOneStationTight',     kstTrkpTMOneStationTight      ,               'kstTrkpTMOneStationTight/F')
ntuple.Branch('kstTrkpTMOneStationLoose',     kstTrkpTMOneStationLoose      ,               'kstTrkpTMOneStationLoose/F')
ntuple.Branch('kstTrkpTrackerMuonArbitrated', kstTrkpTrackerMuonArbitrated  ,               'kstTrkpTrackerMuonArbitrated/F')


ntuple.Branch('bPt',        bPt          ,               'bPt/F')
ntuple.Branch('kstPt',      kstPt        ,               'kstPt/F')
ntuple.Branch('mumuPt',     mumuPt       ,               'mumuPt/F')
ntuple.Branch('mumPt',      mumPt        ,               'mumPt/F')
ntuple.Branch('mupPt',      mupPt        ,               'mupPt/F')
ntuple.Branch('kstTrkmPt',  kstTrkmPt    ,               'kstTrkmPt/F')
ntuple.Branch('kstTrkpPt',  kstTrkpPt    ,               'kstTrkpPt/F')
ntuple.Branch('bPhi',       bPhi         ,               'bPhi/F')
ntuple.Branch('kstPhi',     kstPhi       ,               'kstPhi/F')
ntuple.Branch('mumuPhi',    mumuPhi      ,               'mumuPhi/F')
ntuple.Branch('mumPhi',     mumPhi       ,               'mumPhi/F')
ntuple.Branch('mupPhi',     mupPhi       ,               'mupPhi/F')
ntuple.Branch('kstTrkmPhi', kstTrkmPhi   ,               'kstTrkmPhi/F')
ntuple.Branch('kstTrkpPhi', kstTrkpPhi   ,               'kstTrkpPhi/F')
ntuple.Branch('bEta',       bEta         ,               'bEta/F')
ntuple.Branch('kstEta',     kstEta       ,               'kstEta/F')
ntuple.Branch('mumuEta',    mumuEta      ,               'mumuEta/F')
ntuple.Branch('mumEta',     mumEta       ,               'mumEta/F')
ntuple.Branch('mupEta',     mupEta       ,               'mupEta/F')
ntuple.Branch('kstTrkmEta', kstTrkmEta   ,               'kstTrkmEta/F')
ntuple.Branch('kstTrkpEta', kstTrkpEta   ,               'kstTrkpEta/F')

# ntuple.Branch('mumIsoN'    ,  mumIsoN    )
# ntuple.Branch('mupIsoN'    ,  mupIsoN    )
# ntuple.Branch('kstTrkmIsoN',  kstTrkmIsoN)
# ntuple.Branch('kstTrkpIsoN',  kstTrkpIsoN)


if not skim:
  ntuple.Branch('bCosAlphaVtx',                 bCosAlphaVtx,                                 'bCosAlphaVtx/F')
  ntuple.Branch('bCosAlphaVtxE',                bCosAlphaVtxE,                                'bCosAlphaVtxE/F')
  ntuple.Branch('priVtxCL',                     priVtxCL,                                     'priVtxCL/F')
  ntuple.Branch('priVtxX',                      priVtxX,                                      'priVtxX/F')
  ntuple.Branch('priVtxY',                      priVtxY,                                      'priVtxY/F')
  ntuple.Branch('priVtxZ',                      priVtxZ,                                      'priVtxZ/F')
#   ntuple.Branch('bLVtx',                        bLVtx,                                        'bLVtx/F')
#   ntuple.Branch('bLVtxE',                       bLVtxE,                                       'bLVtxE/F')
#   ntuple.Branch('bDCAVtx',                      bDCAVtx,                                      'bDCAVtx/F')
#   ntuple.Branch('bDCAVtxE',                     bDCAVtxE,                                     'bDCAVtxE/F')
  ntuple.Branch('mumNPixHits',                  mumNPixHits,                                  'mumNPixHits/F')
  ntuple.Branch('mumNTrkHits',                  mumNTrkHits,                                  'mumNTrkHits/F')
  ntuple.Branch('mupNPixHits',                  mupNPixHits,                                  'mupNPixHits/F')
  ntuple.Branch('mupNTrkHits',                  mupNTrkHits,                                  'mupNTrkHits/F')
#   ntuple.Branch('mumKinkChi2',                  mumKinkChi2,                                  'mumKinkChi2/F')
  ntuple.Branch('mumFracHits',                  mumFracHits,                                  'mumFracHits/F')
#   ntuple.Branch('mupKinkChi2',                  mupKinkChi2,                                  'mupKinkChi2/F')
  ntuple.Branch('mupFracHits',                  mupFracHits,                                  'mupFracHits/F')
#   ntuple.Branch('mumMinIP',                     mumMinIP,                                     'mumMinIP/F')
#   ntuple.Branch('mumMinIPE',                    mumMinIPE,                                    'mumMinIPE/F')
#   ntuple.Branch('mupMinIP',                     mupMinIP,                                     'mupMinIP/F')
#   ntuple.Branch('mupMinIPE',                    mupMinIPE,                                    'mupMinIPE/F')
#   ntuple.Branch('kstTrkmDCAVtx',                kstTrkmDCAVtx,                                'kstTrkmDCAVtx/F')
#   ntuple.Branch('kstTrkmDCAVtxE',               kstTrkmDCAVtxE,                               'kstTrkmDCAVtxE/F')
#   ntuple.Branch('kstTrkpDCAVtx',                kstTrkpDCAVtx,                                'kstTrkpDCAVtx/F')
#   ntuple.Branch('kstTrkpDCAVtxE',               kstTrkpDCAVtxE,                               'kstTrkpDCAVtxE/F')
  ntuple.Branch('mumGlobalMuonPromptTight',     mumGlobalMuonPromptTight      ,               'mumGlobalMuonPromptTight/F')
  ntuple.Branch('mumTMLastStationTight',        mumTMLastStationTight         ,               'mumTMLastStationTight/F')
  ntuple.Branch('mumTMLastStationLoose',        mumTMLastStationLoose         ,               'mumTMLastStationLoose/F')
  ntuple.Branch('mumTM2DCompatibilityTight',    mumTM2DCompatibilityTight     ,               'mumTM2DCompatibilityTight/F')
  ntuple.Branch('mumTM2DCompatibilityLoose',    mumTM2DCompatibilityLoose     ,               'mumTM2DCompatibilityLoose/F')
  ntuple.Branch('mumTMLastStationAngTight',     mumTMLastStationAngTight      ,               'mumTMLastStationAngTight/F')
  ntuple.Branch('mumTMLastStationAngLoose',     mumTMLastStationAngLoose      ,               'mumTMLastStationAngLoose/F')
  ntuple.Branch('mumTMOneStationAngTight',      mumTMOneStationAngTight       ,               'mumTMOneStationAngTight/F')
  ntuple.Branch('mumTMOneStationAngLoose',      mumTMOneStationAngLoose       ,               'mumTMOneStationAngLoose/F')
  ntuple.Branch('mumTrackerMuonArbitrated',     mumTrackerMuonArbitrated      ,               'mumTrackerMuonArbitrated/F')
  ntuple.Branch('mupGlobalMuonPromptTight',     mupGlobalMuonPromptTight      ,               'mupGlobalMuonPromptTight/F')
  ntuple.Branch('mupTMLastStationTight',        mupTMLastStationTight         ,               'mupTMLastStationTight/F')
  ntuple.Branch('mupTMLastStationLoose',        mupTMLastStationLoose         ,               'mupTMLastStationLoose/F')
  ntuple.Branch('mupTM2DCompatibilityTight',    mupTM2DCompatibilityTight     ,               'mupTM2DCompatibilityTight/F')
  ntuple.Branch('mupTM2DCompatibilityLoose',    mupTM2DCompatibilityLoose     ,               'mupTM2DCompatibilityLoose/F')
  ntuple.Branch('mupTMLastStationAngTight',     mupTMLastStationAngTight      ,               'mupTMLastStationAngTight/F')
  ntuple.Branch('mupTMLastStationAngLoose',     mupTMLastStationAngLoose      ,               'mupTMLastStationAngLoose/F')
  ntuple.Branch('mupTMOneStationAngTight',      mupTMOneStationAngTight       ,               'mupTMOneStationAngTight/F')
  ntuple.Branch('mupTMOneStationAngLoose',      mupTMOneStationAngLoose       ,               'mupTMOneStationAngLoose/F')
  ntuple.Branch('mupTrackerMuonArbitrated',     mupTrackerMuonArbitrated      ,               'mupTrackerMuonArbitrated/F')
  ntuple.Branch('kstTrkmTM2DCompatibilityTight',kstTrkmTM2DCompatibilityTight ,               'kstTrkmTM2DCompatibilityTight/F')
  ntuple.Branch('kstTrkmTM2DCompatibilityLoose',kstTrkmTM2DCompatibilityLoose ,               'kstTrkmTM2DCompatibilityLoose/F')
  ntuple.Branch('kstTrkpTM2DCompatibilityTight',kstTrkpTM2DCompatibilityTight ,               'kstTrkpTM2DCompatibilityTight/F')
  ntuple.Branch('kstTrkpTM2DCompatibilityLoose',kstTrkpTM2DCompatibilityLoose ,               'kstTrkpTM2DCompatibilityLoose/F')
  ntuple.Branch('kstTrkmTMLastStationAngTight', kstTrkmTMLastStationAngTight  ,               'kstTrkmTMLastStationAngTight/F')
  ntuple.Branch('kstTrkmTMLastStationAngLoose', kstTrkmTMLastStationAngLoose  ,               'kstTrkmTMLastStationAngLoose/F')
  ntuple.Branch('kstTrkpTMOneStationAngTight',  kstTrkpTMOneStationAngTight   ,               'kstTrkpTMOneStationAngTight/F')
  ntuple.Branch('kstTrkpTMOneStationAngLoose',  kstTrkpTMOneStationAngLoose   ,               'kstTrkpTMOneStationAngLoose/F')
  ntuple.Branch('kstTrkmTMOneStationAngTight',  kstTrkmTMOneStationAngTight   ,               'kstTrkmTMOneStationAngTight/F')
  ntuple.Branch('kstTrkmTMOneStationAngLoose',  kstTrkmTMOneStationAngLoose   ,               'kstTrkmTMOneStationAngLoose/F')
  ntuple.Branch('kstTrkpTMLastStationAngTight', kstTrkpTMLastStationAngTight  ,               'kstTrkpTMLastStationAngTight/F')
  ntuple.Branch('kstTrkpTMLastStationAngLoose', kstTrkpTMLastStationAngLoose  ,               'kstTrkpTMLastStationAngLoose/F')
  ntuple.Branch('kstTrkmTMLastStationTight',    kstTrkmTMLastStationTight     ,               'kstTrkmTMLastStationTight/F')
  ntuple.Branch('kstTrkmTMLastStationLoose',    kstTrkmTMLastStationLoose     ,               'kstTrkmTMLastStationLoose/F')
  ntuple.Branch('kstTrkpTMLastStationTight',    kstTrkpTMLastStationTight     ,               'kstTrkpTMLastStationTight/F')
  ntuple.Branch('kstTrkpTMLastStationLoose',    kstTrkpTMLastStationLoose     ,               'kstTrkpTMLastStationLoose/F')
#   ntuple.Branch('kstTrkpMinIP',                 kstTrkpMinIP,                                 'kstTrkpMinIP/F')
#   ntuple.Branch('kstTrkpMinIPE',                kstTrkpMinIPE,                                'kstTrkpMinIPE/F')
#   ntuple.Branch('kstTrkmMinIP',                 kstTrkmMinIP,                                 'kstTrkmMinIP/F')
#   ntuple.Branch('kstTrkmMinIPE',                kstTrkmMinIPE,                                'kstTrkmMinIPE/F')

#   ntuple.Branch('mumIso'    ,  mumIso    )
#   ntuple.Branch('mupIso'    ,  mupIso    )
#   ntuple.Branch('kstTrkmIso',  kstTrkmIso)
#   ntuple.Branch('kstTrkpIso',  kstTrkpIso)


numEvents = tree_lmnr.GetEntries()
print 'total number of events in tree:', numEvents

ROOT.gROOT.LoadMacro('FindValueFromVectorOfBool.h+')

progressbarWidth = 40
sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
sys.stdout.flush()                          # this forces to print the stdout buffer
sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['


for i, ev in enumerate(tree_lmnr):

    if i%int(numEvents/(progressbarWidth-1))==0:
        sys.stdout.write('+')
        sys.stdout.flush()

#     if i > 2000: break

    if not len(ev.bMass) > 0:
        continue
        
    ## trigger requirements
#     if not any( path in ev.TrigTable[0] for path in paths):
#         continue     

    ## now loop on candidates per event
    for icand in range(len(ev.bMass)):
    
        for var in newbr:
            var[0] = -99.
    
        ## muon trigger match: both muons should be matched + 
        ## track trigger match: at least one track should be matched
#         if not ( any( path in ev.mumTrig[icand] and path in ev.mupTrig[icand] and\
#                      (path in ev.kstTrkmTrig[icand] or path in ev.kstTrkpTrig[icand]) for path in paths
#                     ) ):
#             continue

        if skimSoftMu and not (ev.mumNTrkLayers[icand] >= 6        and ev.mupNTrkLayers[icand] >= 6  and \
                               ev.mumNPixLayers[icand] >= 1        and ev.mupNPixLayers[icand] >= 1  and \
#                                ev.mumdxyVtx[icand] < 0.3           and ev.mupdxyVtx[icand] < 0.3     and \
#                                ev.mumdzVtx[icand] < 20             and ev.mupdzVtx[icand]  < 20      and \
                               ROOT.FindValueFromVectorOfBool(ev.mumHighPurity, icand) == 1          and \
                               ROOT.FindValueFromVectorOfBool(ev.mupHighPurity, icand) == 1        ):
            continue


# '''
# special = ['trig']
# for i in branches:
#     if i not in special:
#         print "{NAME}[0]\t = ev.{NAME}[icand]".format(NAME=i).expandtabs(30)
#     if i == 'trig':
#         print "trig[0]                        = paths.index(ev.TrigTable[0].split('_v')[0])"
#     if 'HighPurity' in i:        
#         print "{NAME}[0]\t = ROOT.FindValueFromVectorOfBool(ev.{NAME}, icand)".format(NAME=i).expandtabs(30)
# '''
        ## per event quantities
        runN[0]                        = ev.runN
        eventN[0]                      = ev.eventN
        recoVtxN[0]                    = ev.recoVtxN
        evWeight[0]                    = ev.evWeight
        evWeightE2[0]                  = ev.evWeightE2
        for index,ibx in enumerate(ev.bunchXingMC):
          if ibx==0:
            trueNumInteractionsMC[0]       = ev.trueNumInteractionsMC[index]
            break
#         trig[0]                        = paths.index(ev.TrigTable[0].split('_v')[0])
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
        bPx[0]                         = ev.bPx[icand]
        bPy[0]                         = ev.bPy[icand]
        bPz[0]                         = ev.bPz[icand]
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
#         bctauPVBS[0]                   = ev.bctauPVBS[icand]
#         bctauPVBSE[0]                  = ev.bctauPVBSE[icand]
        kstMass[0]                     = ev.kstMass[icand]
        kstMassE[0]                    = ev.kstMassE[icand]
        kstBarMass[0]                  = ev.kstBarMass[icand]
        kstBarMassE[0]                 = ev.kstBarMassE[icand]
        kstPx[0]                       = ev.kstPx[icand]
        kstPy[0]                       = ev.kstPy[icand]
        kstPz[0]                       = ev.kstPz[icand]
#         kstPxxE[0]                     = ev.kstPxxE[icand]
#         kstPyyE[0]                     = ev.kstPyyE[icand]
#         kstPzzE[0]                     = ev.kstPzzE[icand]
#         kstPxyE[0]                     = ev.kstPxyE[icand]
#         kstPxzE[0]                     = ev.kstPxzE[icand]
#         kstPyzE[0]                     = ev.kstPyzE[icand]
        kstVtxCL[0]                    = ev.kstVtxCL[icand]
        kstVtxX[0]                     = ev.kstVtxX[icand]
        kstVtxY[0]                     = ev.kstVtxY[icand]
        kstVtxZ[0]                     = ev.kstVtxZ[icand]
        kkMass[0]                      = computeInvMass(ev.kstTrkmPx[icand], ev.kstTrkmPy[icand], ev.kstTrkmPz[icand], kaonMass,
                                                        ev.kstTrkpPx[icand], ev.kstTrkpPy[icand], ev.kstTrkpPz[icand], kaonMass )

        
        mumuMass[0]                    = ev.mumuMass[icand]
        mumuMassE[0]                   = ev.mumuMassE[icand]
        mumuPx[0]                      = ev.mumuPx[icand]
        mumuPy[0]                      = ev.mumuPy[icand]
        mumuPz[0]                      = ev.mumuPz[icand]
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
        mumPx[0]                       = ev.mumPx[icand]
        mumPy[0]                       = ev.mumPy[icand]
        mumPz[0]                       = ev.mumPz[icand]
#         mumDCAVtx[0]                   = ev.mumDCAVtx[icand]
#         mumDCAVtxE[0]                  = ev.mumDCAVtxE[icand]
        mumDCABS[0]                    = ev.mumDCABS[icand]
        mumDCABSE[0]                   = ev.mumDCABSE[icand]
#         mumdxyVtx[0]                   = ev.mumdxyVtx[icand]
#         mumdzVtx[0]                    = ev.mumdzVtx[icand]
#         mumMinIP2D[0]                  = ev.mumMinIP2D[icand]
#         mumMinIP2DE[0]                 = ev.mumMinIP2DE[icand]
        mumNPixLayers[0]               = ev.mumNPixLayers[icand]
        mumNTrkLayers[0]               = ev.mumNTrkLayers[icand]
#         mumNMuonHits[0]                = ev.mumNMuonHits[icand]
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
        mupPx[0]                       = ev.mupPx[icand]
        mupPy[0]                       = ev.mupPy[icand]
        mupPz[0]                       = ev.mupPz[icand]
#         mupDCAVtx[0]                   = ev.mupDCAVtx[icand]
#         mupDCAVtxE[0]                  = ev.mupDCAVtxE[icand]
        mupDCABS[0]                    = ev.mupDCABS[icand]
        mupDCABSE[0]                   = ev.mupDCABSE[icand]
#         mupdxyVtx[0]                   = ev.mupdxyVtx[icand]
#         mupdzVtx[0]                    = ev.mupdzVtx[icand]
#         mupMinIP2D[0]                  = ev.mupMinIP2D[icand]
#         mupMinIP2DE[0]                 = ev.mupMinIP2DE[icand]
        mupNPixLayers[0]               = ev.mupNPixLayers[icand]
        mupNTrkLayers[0]               = ev.mupNTrkLayers[icand]
#         mupNMuonHits[0]                = ev.mupNMuonHits[icand]
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
        kstTrkmPx[0]                   = ev.kstTrkmPx[icand]
        kstTrkmPy[0]                   = ev.kstTrkmPy[icand]
        kstTrkmPz[0]                   = ev.kstTrkmPz[icand]
#         kstTrkmPxxE[0]                 = ev.kstTrkmPxxE[icand]
#         kstTrkmPyyE[0]                 = ev.kstTrkmPyyE[icand]
#         kstTrkmPzzE[0]                 = ev.kstTrkmPzzE[icand]
#         kstTrkmPxyE[0]                 = ev.kstTrkmPxyE[icand]
#         kstTrkmPxzE[0]                 = ev.kstTrkmPxzE[icand]
#         kstTrkmPyzE[0]                 = ev.kstTrkmPyzE[icand]
        kstTrkmDCABS[0]                = ev.kstTrkmDCABS[icand]     ## theDCAXBS.perigeeParameters().transverseImpactParameter();
        kstTrkmDCABSE[0]               = ev.kstTrkmDCABSE[icand]
        kstTrkmFracHits[0]             = ev.kstTrkmFracHits[icand]
#         kstTrkmdxyVtx[0]               = ev.kstTrkmdxyVtx[icand]    ## IPTools::absoluteTransverseImpactParameter(TrackmTT, bestVtxReFit); 
#         kstTrkmdzVtx[0]                = ev.kstTrkmdzVtx[icand]     ## TrackmTT.track().dz(bestVtxReFit.position())
        kstTrkmNPixHits[0]             = ev.kstTrkmNPixHits[icand]
        kstTrkmNPixLayers[0]           = ev.kstTrkmNPixLayers[icand]
        kstTrkmNTrkHits[0]             = ev.kstTrkmNTrkHits[icand]
        kstTrkmNTrkLayers[0]           = ev.kstTrkmNTrkLayers[icand]
#         kstTrkmMinIP2D[0]              = ev.kstTrkmMinIP2D[icand]
#         kstTrkmMinIP2DE[0]             = ev.kstTrkmMinIP2DE[icand]

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
        kstTrkpPx[0]                   = ev.kstTrkpPx[icand]
        kstTrkpPy[0]                   = ev.kstTrkpPy[icand]
        kstTrkpPz[0]                   = ev.kstTrkpPz[icand]
#         kstTrkpPxxE[0]                 = ev.kstTrkpPxxE[icand]
#         kstTrkpPyyE[0]                 = ev.kstTrkpPyyE[icand]
#         kstTrkpPzzE[0]                 = ev.kstTrkpPzzE[icand]
#         kstTrkpPxyE[0]                 = ev.kstTrkpPxyE[icand]
#         kstTrkpPxzE[0]                 = ev.kstTrkpPxzE[icand]
#         kstTrkpPyzE[0]                 = ev.kstTrkpPyzE[icand]
        kstTrkpDCABS[0]                = ev.kstTrkpDCABS[icand]
        kstTrkpDCABSE[0]               = ev.kstTrkpDCABSE[icand]
        kstTrkpFracHits[0]             = ev.kstTrkpFracHits[icand]
#         kstTrkpdxyVtx[0]               = ev.kstTrkpdxyVtx[icand]
#         kstTrkpdzVtx[0]                = ev.kstTrkpdzVtx[icand]
        kstTrkpNPixHits[0]             = ev.kstTrkpNPixHits[icand]
        kstTrkpNPixLayers[0]           = ev.kstTrkpNPixLayers[icand]
        kstTrkpNTrkHits[0]             = ev.kstTrkpNTrkHits[icand]
        kstTrkpNTrkLayers[0]           = ev.kstTrkpNTrkLayers[icand]
#         kstTrkpMinIP2D[0]              = ev.kstTrkpMinIP2D[icand]
#         kstTrkpMinIP2DE[0]             = ev.kstTrkpMinIP2DE[icand]

        trkpCategoryDict = muonCategory(ev.kstTrkpMuMatch[icand])
        kstTrkpGlobalMuon            [0] = trkpCategoryDict['GlobalMuon']
        kstTrkpTrackerMuon           [0] = trkpCategoryDict['TrackerMuon']
        kstTrkpStandAloneMuon        [0] = trkpCategoryDict['StandAloneMuon']
        kstTrkpTrackerMuonArbitrated [0] = trkpCategoryDict['TrackerMuonArbitrated']
        kstTrkpTMOneStationTight     [0] = trkpCategoryDict['TMOneStationTight']
        kstTrkpTMOneStationLoose     [0] = trkpCategoryDict['TMOneStationLoose']

        tagB0[0]                      = FlavorTagger(ev.kstMass[icand], ev.kstBarMass[icand])  ## 1 if B0, 0 if B0bar

        
#         mumIsoN    .clear()
#         mupIsoN    .clear()
#         kstTrkmIsoN.clear()
#         kstTrkpIsoN.clear()

#         num = 10
#         for i in range(num):
#             if len( ev.mumIso[icand] ) > 0:
#                 mumIsoN.push_back( numpy.size(numpy.where(numpy.array(ev.mumIso[icand]) < 0.1*(i+1))) ) 
#             else:
#                 mumIsoN.push_back(0) 
# 
#             if len( ev.mupIso[icand] ) > 0:
#                 mupIsoN.push_back( numpy.size(numpy.where(numpy.array(ev.mupIso[icand]) < 0.1*(i+1))) ) 
#             else:
#                 mupIsoN.push_back(0) 
# 
#             if len( ev.kstTrkmIso[icand] ) > 0:
#                 kstTrkmIsoN.push_back( numpy.size(numpy.where(numpy.array(ev.kstTrkmIso[icand]) < 0.1*(i+1))) ) 
#             else:
#                 kstTrkmIsoN.push_back(0) 
# 
#             if len( ev.kstTrkpIso[icand] ) > 0:
#                 kstTrkpIsoN.push_back( numpy.size(numpy.where(numpy.array(ev.kstTrkpIso[icand]) < 0.1*(i+1))) ) 
#             else:
#                 kstTrkpIsoN.push_back(0) 

        if not skim:
#             bCosAlphaVtx[0]                = ev.bCosAlphaVtx[icand]
#             bCosAlphaVtxE[0]               = ev.bCosAlphaVtxE[icand]
            priVtxCL[0]                    = ev.priVtxCL
            priVtxX[0]                     = ev.priVtxX
            priVtxY[0]                     = ev.priVtxY
            priVtxZ[0]                     = ev.priVtxZ
#             bLVtx[0]                       = ev.bLVtx[icand]
#             bLVtxE[0]                      = ev.bLVtxE[icand]
#             bDCAVtx[0]                     = ev.bDCAVtx[icand]
#             bDCAVtxE[0]                    = ev.bDCAVtxE[icand]
            mumNPixHits[0]                 = ev.mumNPixHits[icand]
            mumNTrkHits[0]                 = ev.mumNTrkHits[icand]
#             mumKinkChi2[0]                 = ev.mumKinkChi2[icand]
            mumFracHits[0]                 = ev.mumFracHits[icand]
#             mumMinIP[0]                    = ev.mumMinIP[icand]
#             mumMinIPE[0]                   = ev.mumMinIPE[icand]
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
#             mupKinkChi2[0]                 = ev.mupKinkChi2[icand]
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
#             mupMinIP[0]                    = ev.mupMinIP[icand]
#             mupMinIPE[0]                   = ev.mupMinIPE[icand]
            
#             kstTrkmDCAVtx[0]               = ev.kstTrkmDCAVtx[icand]    ## absoluteImpactParameter3D (track, best PV ReFit )
#             kstTrkmDCAVtxE[0]              = ev.kstTrkmDCAVtxE[icand]
#             kstTrkpDCAVtx[0]               = ev.kstTrkpDCAVtx[icand]
#             kstTrkpDCAVtxE[0]              = ev.kstTrkpDCAVtxE[icand]
            kstTrkmGlobalMuonPromptTight [0] = trkmCategoryDict['GlobalMuonPromptTight']
            kstTrkmTMLastStationTight    [0] = trkmCategoryDict['TMLastStationTight']
            kstTrkmTMLastStationLoose    [0] = trkmCategoryDict['TMLastStationLoose']
            kstTrkmTM2DCompatibilityTight[0] = trkmCategoryDict['TM2DCompatibilityTight']
            kstTrkmTM2DCompatibilityLoose[0] = trkmCategoryDict['TM2DCompatibilityLoose']
            kstTrkmTMLastStationAngTight [0] = trkmCategoryDict['TMLastStationAngTight']
            kstTrkmTMLastStationAngLoose [0] = trkmCategoryDict['TMLastStationAngLoose']
            kstTrkmTMOneStationAngTight  [0] = trkmCategoryDict['TMOneStationAngTight']
            kstTrkmTMOneStationAngLoose  [0] = trkmCategoryDict['TMOneStationAngLoose']
  
#             kstTrkmMinIP[0]                = ev.kstTrkmMinIP[icand]
#             kstTrkmMinIPE[0]               = ev.kstTrkmMinIPE[icand]
#             kstTrkpMinIP[0]                = ev.kstTrkpMinIP[icand]
#             kstTrkpMinIPE[0]               = ev.kstTrkpMinIPE[icand]
            kstTrkpGlobalMuonPromptTight [0] = trkpCategoryDict['GlobalMuonPromptTight']
            kstTrkpTMLastStationTight    [0] = trkpCategoryDict['TMLastStationTight']
            kstTrkpTMLastStationLoose    [0] = trkpCategoryDict['TMLastStationLoose']
            kstTrkpTM2DCompatibilityTight[0] = trkpCategoryDict['TM2DCompatibilityTight']
            kstTrkpTM2DCompatibilityLoose[0] = trkpCategoryDict['TM2DCompatibilityLoose']
            kstTrkpTMLastStationAngTight [0] = trkpCategoryDict['TMLastStationAngTight']
            kstTrkpTMLastStationAngLoose [0] = trkpCategoryDict['TMLastStationAngLoose']
            kstTrkpTMOneStationAngTight  [0] = trkpCategoryDict['TMOneStationAngTight']
            kstTrkpTMOneStationAngLoose  [0] = trkpCategoryDict['TMOneStationAngLoose']

#             mumIso    .clear()
#             mupIso    .clear()
#             kstTrkmIso.clear()
#             kstTrkpIso.clear()
#     
#             if len( ev.mumIso[icand] ) > 0:
#                 for i,isoi in enumerate(ev.mumIso[icand]):
#                     mumIso.push_back(isoi)
#             else:
#                 mumIso.push_back(0)        
#     
#             if len( ev.mupIso[icand] ) > 0:
#                 for i,isoi in enumerate(ev.mupIso[icand]):
#                     mupIso.push_back(isoi)
#             else:
#                 mupIso.push_back(0)        
#     
#             if len( ev.kstTrkmIso[icand] ) > 0:
#                 for i,isoi in enumerate(ev.kstTrkmIso[icand]):
#                     kstTrkmIso.push_back(isoi)
#             else:
#                 kstTrkmIso.push_back(0)        
#     
#             if len( ev.kstTrkpIso[icand] ) > 0:
#                 for i,isoi in enumerate(ev.kstTrkpIso[icand]):
#                     kstTrkpIso.push_back(isoi)
#             else:
#                 kstTrkpIso.push_back(0)        



        if isMC:
            truthMatchMum[0]               = ROOT.FindValueFromVectorOfBool(ev.truthMatchMum,  icand)
            truthMatchMup[0]               = ROOT.FindValueFromVectorOfBool(ev.truthMatchMup,  icand)
            truthMatchTrkm[0]              = ROOT.FindValueFromVectorOfBool(ev.truthMatchTrkm, icand)
            truthMatchTrkp[0]              = ROOT.FindValueFromVectorOfBool(ev.truthMatchTrkp, icand)
            mumDeltaRwithMC[0]             = ev.mumDeltaRwithMC[icand]
            mupDeltaRwithMC[0]             = ev.mupDeltaRwithMC[icand]
            kstTrkpDeltaRwithMC[0]         = ev.kstTrkpDeltaRwithMC[icand]
            kstTrkmDeltaRwithMC[0]         = ev.kstTrkmDeltaRwithMC[icand]

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
            

        ntuple. Fill()


sys.stdout.write('\n')

file_out.cd()
ntuple.Write()
file_out.Close()


