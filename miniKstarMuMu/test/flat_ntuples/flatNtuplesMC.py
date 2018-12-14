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
#  	tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/phaseII/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_PhaseII_KStarMuMu_firstRound_may30/180530_122915/0000/B0ToKstMuMu_*.root')
 	tree_lmnr.Add('/gwteras/cms/store/user/fiorendi/p5prime/phaseII/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_PhaseII_KStarMuMu_firstRound_PU0/180604_075552/0000/B0ToKstMuMu_*.root')

## Mauro 2012, GEN
# tree_lmnr.Add('/gwteras/cms/store/user/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/GENcands/B0ToKstMuMu_GEN_NoFilter_MC_NTuples.root')


file_out = ROOT.TFile('ntuples/flat_ntuple_PhaseII_RECO_%s_PU0.root'%type, 'recreate')
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
    print "{NAME}\t = array('d', [-99.]);  newbr.append({NAME})".format(NAME=i).expandtabs(40)
'''
runN                             = array('d', [-99.]);  newbr.append(runN)
eventN                           = array('d', [-99.]);  newbr.append(eventN)
recoVtxN                         = array('d', [-99.]);  newbr.append(recoVtxN)
evWeight                         = array('d', [-99.]);  newbr.append(evWeight)
evWeightE2                       = array('d', [-99.]);  newbr.append(evWeightE2)
trueNumInteractionsMC            = array('d', [-99.]);  newbr.append(trueNumInteractionsMC)
trig                             = array('d', [-99.]);  newbr.append(trig)
bsX                              = array('d', [-99.]);  newbr.append(bsX)
bsY                              = array('d', [-99.]);  newbr.append(bsY)
bMass                            = array('d', [-99.]);  newbr.append(bMass)
bMassE                           = array('d', [-99.]);  newbr.append(bMassE)
bBarMass                         = array('d', [-99.]);  newbr.append(bBarMass)
bBarMassE                        = array('d', [-99.]);  newbr.append(bBarMassE)
bPx                              = array('d', [-99.]);  newbr.append(bPx)
bPy                              = array('d', [-99.]);  newbr.append(bPy)
bPz                              = array('d', [-99.]);  newbr.append(bPz)
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
# bctauPVBS                        = array('d', [-99.]);  newbr.append(bctauPVBS)
# bctauPVBSE                       = array('d', [-99.]);  newbr.append(bctauPVBSE)
kstMass                          = array('d', [-99.]);  newbr.append(kstMass)
kstMassE                         = array('d', [-99.]);  newbr.append(kstMassE)
kstBarMass                       = array('d', [-99.]);  newbr.append(kstBarMass)
kstBarMassE                      = array('d', [-99.]);  newbr.append(kstBarMassE)
kstPx                            = array('d', [-99.]);  newbr.append(kstPx)
kstPy                            = array('d', [-99.]);  newbr.append(kstPy)
kstPz                            = array('d', [-99.]);  newbr.append(kstPz)
# kstPxxE                          = array('d', [-99.]);  newbr.append(kstPxxE)
# kstPyyE                          = array('d', [-99.]);  newbr.append(kstPyyE)
# kstPzzE                          = array('d', [-99.]);  newbr.append(kstPzzE)
# kstPxyE                          = array('d', [-99.]);  newbr.append(kstPxyE)
# kstPxzE                          = array('d', [-99.]);  newbr.append(kstPxzE)
# kstPyzE                          = array('d', [-99.]);  newbr.append(kstPyzE)
kstVtxCL                         = array('d', [-99.]);  newbr.append(kstVtxCL)
kstVtxX                          = array('d', [-99.]);  newbr.append(kstVtxX)
kstVtxY                          = array('d', [-99.]);  newbr.append(kstVtxY)
kstVtxZ                          = array('d', [-99.]);  newbr.append(kstVtxZ)
kkMass                           = array('d', [-99.]);  newbr.append(kkMass)
mumuMass                         = array('d', [-99.]);  newbr.append(mumuMass)
mumuMassE                        = array('d', [-99.]);  newbr.append(mumuMassE)
mumuPx                           = array('d', [-99.]);  newbr.append(mumuPx)
mumuPy                           = array('d', [-99.]);  newbr.append(mumuPy)
mumuPz                           = array('d', [-99.]);  newbr.append(mumuPz)
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
mumPx                            = array('d', [-99.]);  newbr.append(mumPx)
mumPy                            = array('d', [-99.]);  newbr.append(mumPy)
mumPz                            = array('d', [-99.]);  newbr.append(mumPz)
# mumDCAVtx                        = array('d', [-99.]);  newbr.append(mumDCAVtx)
# mumDCAVtxE                       = array('d', [-99.]);  newbr.append(mumDCAVtxE)
mumDCABS                         = array('d', [-99.]);  newbr.append(mumDCABS)
mumDCABSE                        = array('d', [-99.]);  newbr.append(mumDCABSE)
# mumdxyVtx                        = array('d', [-99.]);  newbr.append(mumdxyVtx)
# mumdzVtx                         = array('d', [-99.]);  newbr.append(mumdzVtx)
# mumMinIP2D                       = array('d', [-99.]);  newbr.append(mumMinIP2D)
# mumMinIP2DE                      = array('d', [-99.]);  newbr.append(mumMinIP2DE)
mumNPixLayers                    = array('d', [-99.]);  newbr.append(mumNPixLayers)
mumNTrkLayers                    = array('d', [-99.]);  newbr.append(mumNTrkLayers)
# mumNMuonHits                     = array('d', [-99.]);  newbr.append(mumNMuonHits)
mumNMatchStation                 = array('d', [-99.]);  newbr.append(mumNMatchStation)
mupHighPurity                    = array('d', [-99.]);  newbr.append(mupHighPurity)
mupCL                            = array('d', [-99.]);  newbr.append(mupCL)
mupNormChi2                      = array('d', [-99.]);  newbr.append(mupNormChi2)
mupPx                            = array('d', [-99.]);  newbr.append(mupPx)
mupPy                            = array('d', [-99.]);  newbr.append(mupPy)
mupPz                            = array('d', [-99.]);  newbr.append(mupPz)
# mupDCAVtx                        = array('d', [-99.]);  newbr.append(mupDCAVtx)
# mupDCAVtxE                       = array('d', [-99.]);  newbr.append(mupDCAVtxE)
mupDCABS                         = array('d', [-99.]);  newbr.append(mupDCABS)
mupDCABSE                        = array('d', [-99.]);  newbr.append(mupDCABSE)
mupdxyVtx                        = array('d', [-99.]);  newbr.append(mupdxyVtx)
mupdzVtx                         = array('d', [-99.]);  newbr.append(mupdzVtx)
# mupMinIP2D                       = array('d', [-99.]);  newbr.append(mupMinIP2D)
# mupMinIP2DE                      = array('d', [-99.]);  newbr.append(mupMinIP2DE)
mupNPixLayers                    = array('d', [-99.]);  newbr.append(mupNPixLayers)
mupNTrkLayers                    = array('d', [-99.]);  newbr.append(mupNTrkLayers)
# mupNMuonHits                     = array('d', [-99.]);  newbr.append(mupNMuonHits)
mupNMatchStation                 = array('d', [-99.]);  newbr.append(mupNMatchStation)
kstTrkmHighPurity                = array('d', [-99.]);  newbr.append(kstTrkmHighPurity)
kstTrkmCL                        = array('d', [-99.]);  newbr.append(kstTrkmCL)
kstTrkmNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkmNormChi2)
kstTrkmPx                        = array('d', [-99.]);  newbr.append(kstTrkmPx)
kstTrkmPy                        = array('d', [-99.]);  newbr.append(kstTrkmPy)
kstTrkmPz                        = array('d', [-99.]);  newbr.append(kstTrkmPz)
# kstTrkmPxxE                      = array('d', [-99.]);  newbr.append(kstTrkmPxxE)
# kstTrkmPyyE                      = array('d', [-99.]);  newbr.append(kstTrkmPyyE)
# kstTrkmPzzE                      = array('d', [-99.]);  newbr.append(kstTrkmPzzE)
# kstTrkmPxyE                      = array('d', [-99.]);  newbr.append(kstTrkmPxyE)
# kstTrkmPxzE                      = array('d', [-99.]);  newbr.append(kstTrkmPxzE)
# kstTrkmPyzE                      = array('d', [-99.]);  newbr.append(kstTrkmPyzE)
kstTrkmDCABS                     = array('d', [-99.]);  newbr.append(kstTrkmDCABS)
kstTrkmDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkmDCABSE)
kstTrkmFracHits                  = array('d', [-99.]);  newbr.append(kstTrkmFracHits)
kstTrkmdxyVtx                    = array('d', [-99.]);  newbr.append(kstTrkmdxyVtx)
kstTrkmdzVtx                     = array('d', [-99.]);  newbr.append(kstTrkmdzVtx)
# kstTrkmMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkmMinIP2D)
# kstTrkmMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkmMinIP2DE)
kstTrkmNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkmNPixHits)
kstTrkmNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkmNPixLayers)
kstTrkmNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkmNTrkHits)
kstTrkmNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkmNTrkLayers)
kstTrkmMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkmMuMatch)
kstTrkpHighPurity                = array('d', [-99.]);  newbr.append(kstTrkpHighPurity)
kstTrkpCL                        = array('d', [-99.]);  newbr.append(kstTrkpCL)
kstTrkpNormChi2                  = array('d', [-99.]);  newbr.append(kstTrkpNormChi2)
kstTrkpPx                        = array('d', [-99.]);  newbr.append(kstTrkpPx)
kstTrkpPy                        = array('d', [-99.]);  newbr.append(kstTrkpPy)
kstTrkpPz                        = array('d', [-99.]);  newbr.append(kstTrkpPz)
# kstTrkpPxxE                      = array('d', [-99.]);  newbr.append(kstTrkpPxxE)
# kstTrkpPyyE                      = array('d', [-99.]);  newbr.append(kstTrkpPyyE)
# kstTrkpPzzE                      = array('d', [-99.]);  newbr.append(kstTrkpPzzE)
# kstTrkpPxyE                      = array('d', [-99.]);  newbr.append(kstTrkpPxyE)
# kstTrkpPxzE                      = array('d', [-99.]);  newbr.append(kstTrkpPxzE)
# kstTrkpPyzE                      = array('d', [-99.]);  newbr.append(kstTrkpPyzE)
kstTrkpDCABS                     = array('d', [-99.]);  newbr.append(kstTrkpDCABS)
kstTrkpDCABSE                    = array('d', [-99.]);  newbr.append(kstTrkpDCABSE)
kstTrkpFracHits                  = array('d', [-99.]);  newbr.append(kstTrkpFracHits)
kstTrkpdxyVtx                    = array('d', [-99.]);  newbr.append(kstTrkpdxyVtx)
kstTrkpdzVtx                     = array('d', [-99.]);  newbr.append(kstTrkpdzVtx)
# kstTrkpMinIP2D                   = array('d', [-99.]);  newbr.append(kstTrkpMinIP2D)
# kstTrkpMinIP2DE                  = array('d', [-99.]);  newbr.append(kstTrkpMinIP2DE)
kstTrkpNPixHits                  = array('d', [-99.]);  newbr.append(kstTrkpNPixHits)
kstTrkpNPixLayers                = array('d', [-99.]);  newbr.append(kstTrkpNPixLayers)
kstTrkpNTrkHits                  = array('d', [-99.]);  newbr.append(kstTrkpNTrkHits)
kstTrkpNTrkLayers                = array('d', [-99.]);  newbr.append(kstTrkpNTrkLayers)
kstTrkpMuMatch                   = array('d', [-99.]);  newbr.append(kstTrkpMuMatch)
tagB0                            = array('d', [-99.]);  newbr.append(tagB0)

truthMatchMum                    = array('d', [-99.]);  newbr.append(truthMatchMum)
truthMatchMup                    = array('d', [-99.]);  newbr.append(truthMatchMup)
truthMatchTrkm                   = array('d', [-99.]);  newbr.append(truthMatchTrkm)
truthMatchTrkp                   = array('d', [-99.]);  newbr.append(truthMatchTrkp)
mumDeltaRwithMC                  = array('d', [-99.]);  newbr.append(mumDeltaRwithMC)
mupDeltaRwithMC                  = array('d', [-99.]);  newbr.append(mupDeltaRwithMC)
kstTrkpDeltaRwithMC              = array('d', [-99.]);  newbr.append(kstTrkpDeltaRwithMC)
kstTrkmDeltaRwithMC              = array('d', [-99.]);  newbr.append(kstTrkmDeltaRwithMC)
genSignal                        = array('d', [-99.]);  newbr.append(genSignal)
genSignHasFSR                    = array('d', [-99.]);  newbr.append(genSignHasFSR)
genSignKstHasFSR                 = array('d', [-99.]);  newbr.append(genSignKstHasFSR)
genSignPsiHasFSR                 = array('d', [-99.]);  newbr.append(genSignPsiHasFSR)
genPriVtxX                       = array('d', [-99.]);  newbr.append(genPriVtxX)
genPriVtxY                       = array('d', [-99.]);  newbr.append(genPriVtxY)
genPriVtxZ                       = array('d', [-99.]);  newbr.append(genPriVtxZ)
genB0Mass                        = array('d', [-99.]);  newbr.append(genB0Mass)
genB0Px                          = array('d', [-99.]);  newbr.append(genB0Px)
genB0Py                          = array('d', [-99.]);  newbr.append(genB0Py)
genB0Pz                          = array('d', [-99.]);  newbr.append(genB0Pz)
genB0VtxX                        = array('d', [-99.]);  newbr.append(genB0VtxX)
genB0VtxY                        = array('d', [-99.]);  newbr.append(genB0VtxY)
genB0VtxZ                        = array('d', [-99.]);  newbr.append(genB0VtxZ)
genKstMass                       = array('d', [-99.]);  newbr.append(genKstMass)
genKstPx                         = array('d', [-99.]);  newbr.append(genKstPx)
genKstPy                         = array('d', [-99.]);  newbr.append(genKstPy)
genKstPz                         = array('d', [-99.]);  newbr.append(genKstPz)
genKstVtxX                       = array('d', [-99.]);  newbr.append(genKstVtxX)
genKstVtxY                       = array('d', [-99.]);  newbr.append(genKstVtxY)
genKstVtxZ                       = array('d', [-99.]);  newbr.append(genKstVtxZ)
genPsiMass                       = array('d', [-99.]);  newbr.append(genPsiMass)
genPsiVtxX                       = array('d', [-99.]);  newbr.append(genPsiVtxX)
genPsiVtxY                       = array('d', [-99.]);  newbr.append(genPsiVtxY)
genPsiVtxZ                       = array('d', [-99.]);  newbr.append(genPsiVtxZ)
genMumMother                     = array('d', [-99.]);  newbr.append(genMumMother)
genMumPx                         = array('d', [-99.]);  newbr.append(genMumPx)
genMumPy                         = array('d', [-99.]);  newbr.append(genMumPy)
genMumPz                         = array('d', [-99.]);  newbr.append(genMumPz)
genMupMother                     = array('d', [-99.]);  newbr.append(genMupMother)
genMupPx                         = array('d', [-99.]);  newbr.append(genMupPx)
genMupPy                         = array('d', [-99.]);  newbr.append(genMupPy)
genMupPz                         = array('d', [-99.]);  newbr.append(genMupPz)
genKstTrkmMother                 = array('d', [-99.]);  newbr.append(genKstTrkmMother)
genKstTrkmID                     = array('d', [-99.]);  newbr.append(genKstTrkmID)
genKstTrkmPx                     = array('d', [-99.]);  newbr.append(genKstTrkmPx)
genKstTrkmPy                     = array('d', [-99.]);  newbr.append(genKstTrkmPy)
genKstTrkmPz                     = array('d', [-99.]);  newbr.append(genKstTrkmPz)
genKstTrkpMother                 = array('d', [-99.]);  newbr.append(genKstTrkpMother)
genKstTrkpID                     = array('d', [-99.]);  newbr.append(genKstTrkpID)
genKstTrkpPx                     = array('d', [-99.]);  newbr.append(genKstTrkpPx)
genKstTrkpPy                     = array('d', [-99.]);  newbr.append(genKstTrkpPy)
genKstTrkpPz                     = array('d', [-99.]);  newbr.append(genKstTrkpPz)

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


# mumIsoN                          = ROOT.std.vector('float')()   
# mupIsoN                          = ROOT.std.vector('float')()   
# kstTrkmIsoN                      = ROOT.std.vector('float')()   
# kstTrkpIsoN                      = ROOT.std.vector('float')()   

if not skim:
  bCosAlphaVtx                     = array('d', [-99.]);  newbr.append(bCosAlphaVtx)
  bCosAlphaVtxE                    = array('d', [-99.]);  newbr.append(bCosAlphaVtxE)
  priVtxCL                         = array('d', [-99.]);  newbr.append(priVtxCL)
  priVtxX                          = array('d', [-99.]);  newbr.append(priVtxX)
  priVtxY                          = array('d', [-99.]);  newbr.append(priVtxY)
  priVtxZ                          = array('d', [-99.]);  newbr.append(priVtxZ)
#   bLVtx                            = array('d', [-99.]);  newbr.append(bLVtx)
#   bLVtxE                           = array('d', [-99.]);  newbr.append(bLVtxE)
#   bDCAVtx                          = array('d', [-99.]);  newbr.append(bDCAVtx)
#   bDCAVtxE                         = array('d', [-99.]);  newbr.append(bDCAVtxE)
  mumNPixHits                      = array('d', [-99.]);  newbr.append(mumNPixHits)
  mumNTrkHits                      = array('d', [-99.]);  newbr.append(mumNTrkHits)
  mupNPixHits                      = array('d', [-99.]);  newbr.append(mupNPixHits)
  mupNTrkHits                      = array('d', [-99.]);  newbr.append(mupNTrkHits)
#   mumKinkChi2                      = array('d', [-99.]);  newbr.append(mumKinkChi2)
  mumFracHits                      = array('d', [-99.]);  newbr.append(mumFracHits)
#   mupKinkChi2                      = array('d', [-99.]);  newbr.append(mupKinkChi2)
  mupFracHits                      = array('d', [-99.]);  newbr.append(mupFracHits)
#   mumMinIP                         = array('d', [-99.]);  newbr.append(mumMinIP)
#   mumMinIPE                        = array('d', [-99.]);  newbr.append(mumMinIPE)
#   mupMinIP                         = array('d', [-99.]);  newbr.append(mupMinIP)
#   mupMinIPE                        = array('d', [-99.]);  newbr.append(mupMinIPE)
#   kstTrkmDCAVtx                    = array('d', [-99.]);  newbr.append(kstTrkmDCAVtx)
#   kstTrkmDCAVtxE                   = array('d', [-99.]);  newbr.append(kstTrkmDCAVtxE)
#   kstTrkpDCAVtx                    = array('d', [-99.]);  newbr.append(kstTrkpDCAVtx)
#   kstTrkpDCAVtxE                   = array('d', [-99.]);  newbr.append(kstTrkpDCAVtxE)
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
#   kstTrkmMinIP                     = array('d', [-99.]);  newbr.append(kstTrkmMinIP)
#   kstTrkmMinIPE                    = array('d', [-99.]);  newbr.append(kstTrkmMinIPE)
#   kstTrkpMinIP                     = array('d', [-99.]);  newbr.append(kstTrkpMinIP)
#   kstTrkpMinIPE                    = array('d', [-99.]);  newbr.append(kstTrkpMinIPE)

#   mumIso                           = ROOT.std.vector('float')()  
#   mupIso                           = ROOT.std.vector('float')()  
#   kstTrkmIso                       = ROOT.std.vector('float')()  
#   kstTrkpIso                       = ROOT.std.vector('float')()  


'''
for i in branches:
    print "ntuple.Branch('{NAME}',\t{NAME},\t'{NAME}/F')".format(NAME=i).expandtabs(40)
'''
ntuple.Branch('runN',                   runN,                                   'runN/D')
ntuple.Branch('eventN',                 eventN,                                 'eventN/D')
ntuple.Branch('recoVtxN',               recoVtxN,                               'recoVtxN/D')
ntuple.Branch('evWeight',               evWeight,                               'evWeight/D')
ntuple.Branch('evWeightE2',             evWeightE2,                             'evWeightE2/D')
ntuple.Branch('trueNumInteractionsMC',  trueNumInteractionsMC,                  'trueNumInteractionsMC/D')
ntuple.Branch('trig',                   trig,                                   'trig/D')
ntuple.Branch('bsX',                    bsX,                                    'bsX/D')
ntuple.Branch('bsY',                    bsY,                                    'bsY/D')
ntuple.Branch('bMass',                  bMass,                                  'bMass/D')
ntuple.Branch('bMassE',                 bMassE,                                 'bMassE/D')
ntuple.Branch('bBarMass',               bBarMass,                               'bBarMass/D')
ntuple.Branch('bBarMassE',              bBarMassE,                              'bBarMassE/D')
ntuple.Branch('bPx',                    bPx,                                    'bPx/D')
ntuple.Branch('bPy',                    bPy,                                    'bPy/D')
ntuple.Branch('bPz',                    bPz,                                    'bPz/D')
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
# ntuple.Branch('bctauPVBS',              bctauPVBS,                              'bctauPVBS/D')
# ntuple.Branch('bctauPVBSE',             bctauPVBSE,                             'bctauPVBSE/D')
ntuple.Branch('kstMass',                kstMass,                                'kstMass/D')
ntuple.Branch('kstMassE',               kstMassE,                               'kstMassE/D')
ntuple.Branch('kstBarMass',             kstBarMass,                             'kstBarMass/D')
ntuple.Branch('kstBarMassE',            kstBarMassE,                            'kstBarMassE/D')
ntuple.Branch('kstPx',                  kstPx,                                  'kstPx/D')
ntuple.Branch('kstPy',                  kstPy,                                  'kstPy/D')
ntuple.Branch('kstPz',                  kstPz,                                  'kstPz/D')
# ntuple.Branch('kstPxxE',                kstPxxE,                                'kstPxxE/D')
# ntuple.Branch('kstPyyE',                kstPyyE,                                'kstPyyE/D')
# ntuple.Branch('kstPzzE',                kstPzzE,                                'kstPzzE/D')
# ntuple.Branch('kstPxyE',                kstPxyE,                                'kstPxyE/D')
# ntuple.Branch('kstPxzE',                kstPxzE,                                'kstPxzE/D')
# ntuple.Branch('kstPyzE',                kstPyzE,                                'kstPyzE/D')
ntuple.Branch('kstVtxCL',               kstVtxCL,                               'kstVCL/D')
ntuple.Branch('kstVtxX',                kstVtxX,                                'kstVtxX/D')
ntuple.Branch('kstVtxY',                kstVtxY,                                'kstVtxY/D')
ntuple.Branch('kstVtxZ',                kstVtxZ,                                'kstVtxZ/D')
ntuple.Branch('kkMass',                 kkMass,                                 'kkMass/D')
ntuple.Branch('mumuMass',               mumuMass,                               'mumuMass/D')
ntuple.Branch('mumuMassE',              mumuMassE,                              'mumuMassE/D')
ntuple.Branch('mumuPx',                 mumuPx,                                 'mumuPx/D')
ntuple.Branch('mumuPy',                 mumuPy,                                 'mumuPy/D')
ntuple.Branch('mumuPz',                 mumuPz,                                 'mumuPz/D')
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
ntuple.Branch('mumPx',                  mumPx,                                  'mumPx/D')
ntuple.Branch('mumPy',                  mumPy,                                  'mumPy/D')
ntuple.Branch('mumPz',                  mumPz,                                  'mumPz/D')
# ntuple.Branch('mumDCAVtx',              mumDCAVtx,                              'mumDCAVtx/D')
# ntuple.Branch('mumDCAVtxE',             mumDCAVtxE,                             'mumDCAVtxE/D')
ntuple.Branch('mumDCABS',               mumDCABS,                               'mumDCABS/D')
ntuple.Branch('mumDCABSE',              mumDCABSE,                              'mumDCABSE/D')
# ntuple.Branch('mumdxyVtx',              mumdxyVtx,                              'mumdxyVtx/D')
# ntuple.Branch('mumdzVtx',               mumdzVtx,                               'mumdzVtx/D')
# ntuple.Branch('mumMinIP2D',             mumMinIP2D,                             'mumMinIP2D/D')
# ntuple.Branch('mumMinIP2DE',            mumMinIP2DE,                            'mumMinIP2DE/D')
ntuple.Branch('mumNPixLayers',          mumNPixLayers,                          'mumNPixLayers/D')
ntuple.Branch('mumNTrkLayers',          mumNTrkLayers,                          'mumNTrkLayers/D')
# ntuple.Branch('mumNMuonHits',           mumNMuonHits,                           'mumNMuonHits/D')
ntuple.Branch('mumNMatchStation',       mumNMatchStation,                       'mumNMatchStation/D')
ntuple.Branch('mupHighPurity',          mupHighPurity,                          'mupHighPurity/D')
ntuple.Branch('mupCL',                  mupCL,                                  'mupCL/D')
ntuple.Branch('mupNormChi2',            mupNormChi2,                            'mupNormChi2/D')
ntuple.Branch('mupPx',                  mupPx,                                  'mupPx/D')
ntuple.Branch('mupPy',                  mupPy,                                  'mupPy/D')
ntuple.Branch('mupPz',                  mupPz,                                  'mupPz/D')
# ntuple.Branch('mupDCAVtx',              mupDCAVtx,                              'mupDCAVtx/D')
# ntuple.Branch('mupDCAVtxE',             mupDCAVtxE,                             'mupDCAVtxE/D')
ntuple.Branch('mupDCABS',               mupDCABS,                               'mupDCABS/D')
ntuple.Branch('mupDCABSE',              mupDCABSE,                              'mupDCABSE/D')
ntuple.Branch('mupdxyVtx',              mupdxyVtx,                              'mupdxyVtx/D')
ntuple.Branch('mupdzVtx',               mupdzVtx,                               'mupdzVtx/D')
# ntuple.Branch('mupMinIP2D',             mupMinIP2D,                             'mupMinIP2D/D')
# ntuple.Branch('mupMinIP2DE',            mupMinIP2DE,                            'mupMinIP2DE/D')
ntuple.Branch('mupNPixLayers',          mupNPixLayers,                          'mupNPixLayers/D')
ntuple.Branch('mupNTrkLayers',          mupNTrkLayers,                          'mupNTrkLayers/D')
# ntuple.Branch('mupNMuonHits',           mupNMuonHits,                           'mupNMuonHits/D')
ntuple.Branch('mupNMatchStation',       mupNMatchStation,                       'mupNMatchStation/D')

ntuple.Branch('kstTrkmHighPurity',      kstTrkmHighPurity,                      'kstTrkmHighPurity/D')
ntuple.Branch('kstTrkmCL',              kstTrkmCL,                              'kstTrkmCL/D')
ntuple.Branch('kstTrkmNormChi2',        kstTrkmNormChi2,                        'kstTrkmNormChi2/D')
ntuple.Branch('kstTrkmPx',              kstTrkmPx,                              'kstTrkmPx/D')
ntuple.Branch('kstTrkmPy',              kstTrkmPy,                              'kstTrkmPy/D')
ntuple.Branch('kstTrkmPz',              kstTrkmPz,                              'kstTrkmPz/D')
# ntuple.Branch('kstTrkmPxxE',            kstTrkmPxxE,                            'kstTrkmPxxE/D')
# ntuple.Branch('kstTrkmPyyE',            kstTrkmPyyE,                            'kstTrkmPyyE/D')
# ntuple.Branch('kstTrkmPzzE',            kstTrkmPzzE,                            'kstTrkmPzzE/D')
# ntuple.Branch('kstTrkmPxyE',            kstTrkmPxyE,                            'kstTrkmPxyE/D')
# ntuple.Branch('kstTrkmPxzE',            kstTrkmPxzE,                            'kstTrkmPxzE/D')
# ntuple.Branch('kstTrkmPyzE',            kstTrkmPyzE,                            'kstTrkmPyzE/D')
ntuple.Branch('kstTrkmDCABS',           kstTrkmDCABS,                           'kstTrkmDCABS/D')
ntuple.Branch('kstTrkmDCABSE',          kstTrkmDCABSE,                          'kstTrkmDCABSE/D')
ntuple.Branch('kstTrkmFracHits',        kstTrkmFracHits,                        'kstTrkmFracHits/D')
ntuple.Branch('kstTrkmdxyVtx',          kstTrkmdxyVtx,                          'kstTrkmdxyVtx/D')
ntuple.Branch('kstTrkmdzVtx',           kstTrkmdzVtx,                           'kstTrkmdzVtx/D')
# ntuple.Branch('kstTrkmMinIP2D',         kstTrkmMinIP2D,                         'kstTrkmMinIP2D/D')
# ntuple.Branch('kstTrkmMinIP2DE',        kstTrkmMinIP2DE,                        'kstTrkmMinIP2DE/D')
ntuple.Branch('kstTrkmNPixHits',        kstTrkmNPixHits,                        'kstTrkmNPixHits/D')
ntuple.Branch('kstTrkmNPixLayers',      kstTrkmNPixLayers,                      'kstTrkmNPixLayers/D')
ntuple.Branch('kstTrkmNTrkHits',        kstTrkmNTrkHits,                        'kstTrkmNTrkHits/D')
ntuple.Branch('kstTrkmNTrkLayers',      kstTrkmNTrkLayers,                      'kstTrkmNTrkLayers/D')
ntuple.Branch('kstTrkmMuMatch',         kstTrkmMuMatch,                         'kstTrkmMuMatch/D')
ntuple.Branch('kstTrkpHighPurity',      kstTrkpHighPurity,                      'kstTrkpHighPurity/D')
ntuple.Branch('kstTrkpCL',              kstTrkpCL,                              'kstTrkpCL/D')
ntuple.Branch('kstTrkpNormChi2',        kstTrkpNormChi2,                        'kstTrkpNormChi2/D')
ntuple.Branch('kstTrkpPx',              kstTrkpPx,                              'kstTrkpPx/D')
ntuple.Branch('kstTrkpPy',              kstTrkpPy,                              'kstTrkpPy/D')
ntuple.Branch('kstTrkpPz',              kstTrkpPz,                              'kstTrkpPz/D')
# ntuple.Branch('kstTrkpPxxE',            kstTrkpPxxE,                            'kstTrkpPxxE/D')
# ntuple.Branch('kstTrkpPyyE',            kstTrkpPyyE,                            'kstTrkpPyyE/D')
# ntuple.Branch('kstTrkpPzzE',            kstTrkpPzzE,                            'kstTrkpPzzE/D')
# ntuple.Branch('kstTrkpPxyE',            kstTrkpPxyE,                            'kstTrkpPxyE/D')
# ntuple.Branch('kstTrkpPxzE',            kstTrkpPxzE,                            'kstTrkpPxzE/D')
# ntuple.Branch('kstTrkpPyzE',            kstTrkpPyzE,                            'kstTrkpPyzE/D')
ntuple.Branch('kstTrkpDCABS',           kstTrkpDCABS,                           'kstTrkpDCABS/D')
ntuple.Branch('kstTrkpDCABSE',          kstTrkpDCABSE,                          'kstTrkpDCABSE/D')
ntuple.Branch('kstTrkpFracHits',        kstTrkpFracHits,                        'kstTrkpFracHits/D')
ntuple.Branch('kstTrkpdxyVtx',          kstTrkpdxyVtx,                          'kstTrkpdxyVtx/D')
ntuple.Branch('kstTrkpdzVtx',           kstTrkpdzVtx,                           'kstTrkpdzVtx/D')
# ntuple.Branch('kstTrkpMinIP2D',         kstTrkpMinIP2D,                         'kstTrkpMinIP2D/D')
# ntuple.Branch('kstTrkpMinIP2DE',        kstTrkpMinIP2DE,                        'kstTrkpMinIP2DE/D')
ntuple.Branch('kstTrkpNPixHits',        kstTrkpNPixHits,                        'kstTrkpNPixHits/D')
ntuple.Branch('kstTrkpNPixLayers',      kstTrkpNPixLayers,                      'kstTrkpNPixLayers/D')
ntuple.Branch('kstTrkpNTrkHits',        kstTrkpNTrkHits,                        'kstTrkpNTrkHits/D')
ntuple.Branch('kstTrkpNTrkLayers',      kstTrkpNTrkLayers,                      'kstTrkpNTrkLayers/D')
ntuple.Branch('kstTrkpMuMatch',         kstTrkpMuMatch,                         'kstTrkpMuMatch/D')
ntuple.Branch('tagB0',                  tagB0,                                  'tagB0/D')
ntuple.Branch('truthMatchMum',          truthMatchMum,                          'truthMatchMum/D')
ntuple.Branch('truthMatchMup',          truthMatchMup,                          'truthMatchMup/D')
ntuple.Branch('truthMatchTrkm',         truthMatchTrkm,                         'truthMatchTrkm/D')
ntuple.Branch('truthMatchTrkp',         truthMatchTrkp,                         'truthMatchTrkp/D')
ntuple.Branch('mumDeltaRwithMC',        mumDeltaRwithMC,                        'mumDeltaRwithMC/D')
ntuple.Branch('mupDeltaRwithMC',        mupDeltaRwithMC,                        'mupDeltaRwithMC/D')
ntuple.Branch('kstTrkpDeltaRwithMC',    kstTrkpDeltaRwithMC,                    'kstTrkpDeltaRwithMC/D')
ntuple.Branch('kstTrkmDeltaRwithMC',    kstTrkmDeltaRwithMC,                    'kstTrkmDeltaRwithMC/D')
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

# ntuple.Branch('mumIsoN'    ,  mumIsoN    )
# ntuple.Branch('mupIsoN'    ,  mupIsoN    )
# ntuple.Branch('kstTrkmIsoN',  kstTrkmIsoN)
# ntuple.Branch('kstTrkpIsoN',  kstTrkpIsoN)


if not skim:
  ntuple.Branch('bCosAlphaVtx',                 bCosAlphaVtx,                                 'bCosAlphaVtx/D')
  ntuple.Branch('bCosAlphaVtxE',                bCosAlphaVtxE,                                'bCosAlphaVtxE/D')
  ntuple.Branch('priVtxCL',                     priVtxCL,                                     'priVtxCL/D')
  ntuple.Branch('priVtxX',                      priVtxX,                                      'priVtxX/D')
  ntuple.Branch('priVtxY',                      priVtxY,                                      'priVtxY/D')
  ntuple.Branch('priVtxZ',                      priVtxZ,                                      'priVtxZ/D')
#   ntuple.Branch('bLVtx',                        bLVtx,                                        'bLVtx/D')
#   ntuple.Branch('bLVtxE',                       bLVtxE,                                       'bLVtxE/D')
#   ntuple.Branch('bDCAVtx',                      bDCAVtx,                                      'bDCAVtx/D')
#   ntuple.Branch('bDCAVtxE',                     bDCAVtxE,                                     'bDCAVtxE/D')
  ntuple.Branch('mumNPixHits',                  mumNPixHits,                                  'mumNPixHits/D')
  ntuple.Branch('mumNTrkHits',                  mumNTrkHits,                                  'mumNTrkHits/D')
  ntuple.Branch('mupNPixHits',                  mupNPixHits,                                  'mupNPixHits/D')
  ntuple.Branch('mupNTrkHits',                  mupNTrkHits,                                  'mupNTrkHits/D')
#   ntuple.Branch('mumKinkChi2',                  mumKinkChi2,                                  'mumKinkChi2/D')
  ntuple.Branch('mumFracHits',                  mumFracHits,                                  'mumFracHits/D')
#   ntuple.Branch('mupKinkChi2',                  mupKinkChi2,                                  'mupKinkChi2/D')
  ntuple.Branch('mupFracHits',                  mupFracHits,                                  'mupFracHits/D')
#   ntuple.Branch('mumMinIP',                     mumMinIP,                                     'mumMinIP/D')
#   ntuple.Branch('mumMinIPE',                    mumMinIPE,                                    'mumMinIPE/D')
#   ntuple.Branch('mupMinIP',                     mupMinIP,                                     'mupMinIP/D')
#   ntuple.Branch('mupMinIPE',                    mupMinIPE,                                    'mupMinIPE/D')
#   ntuple.Branch('kstTrkmDCAVtx',                kstTrkmDCAVtx,                                'kstTrkmDCAVtx/D')
#   ntuple.Branch('kstTrkmDCAVtxE',               kstTrkmDCAVtxE,                               'kstTrkmDCAVtxE/D')
#   ntuple.Branch('kstTrkpDCAVtx',                kstTrkpDCAVtx,                                'kstTrkpDCAVtx/D')
#   ntuple.Branch('kstTrkpDCAVtxE',               kstTrkpDCAVtxE,                               'kstTrkpDCAVtxE/D')
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
#   ntuple.Branch('kstTrkpMinIP',                 kstTrkpMinIP,                                 'kstTrkpMinIP/D')
#   ntuple.Branch('kstTrkpMinIPE',                kstTrkpMinIPE,                                'kstTrkpMinIPE/D')
#   ntuple.Branch('kstTrkmMinIP',                 kstTrkmMinIP,                                 'kstTrkmMinIP/D')
#   ntuple.Branch('kstTrkmMinIPE',                kstTrkmMinIPE,                                'kstTrkmMinIPE/D')

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


