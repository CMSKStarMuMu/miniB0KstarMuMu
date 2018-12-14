import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("-d", "--doubleg"  , dest = "doubleg", help = "Define if use 1 or 2 gaussian for the signal model", default = 1)
parser.add_argument("-b", "--use_bdt"  , dest = "use_bdt", help = "Define if use BDT (1) or cut n count (0) selection", default = 1)

args = parser.parse_args()


'''
code to fit the B0 mass distribution:
- unbinned fit
- properly take into account the CP state assignment when plotting the mass variable
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
'''

import ROOT
from ROOT import gSystem

gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev
import sys
import math	

B0Mass_   = 5.27958
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

nSigma_psiRej = 3.


tree = ROOT.TChain('ntuple')
for filename in args.inputfile.split(' '):
    print(filename)
    tree.AddFile(filename)


bMass     = RooRealVar("bMass"    , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
bBarMass  = RooRealVar("bBarMass" , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
mumuMass  = RooRealVar("mumuMass" , "mumuMass" , 0, 6);
mumuMassE = RooRealVar("mumuMassE", "mumuMassE", 0, 10000);
tagB0     = RooRealVar("tagB0"    , "tagB0"    , 0, 2);


B0Mass     = RooRealVar("B0Mass"    , "B0Mass"  , B0Mass_  )
JPsiMass   = RooRealVar("JPsiMass"  , "JPsiMass", JPsiMass_)
PsiPMass   = RooRealVar("PsiPMass"  , "PsiPMass", PsiPMass_)
KStMass    = RooRealVar("KStMass"   , "KStMass" , KStMass_ )


thevars = RooArgSet()
thevars.add(bMass)
thevars.add(bBarMass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)

fulldata   = RooDataSet('fulldata', 'fulldataset', tree,  RooArgSet(thevars))

## add to the input tree the combination of the variables for the B0 arb. mass
theBMassfunc = RooFormulaVar("theBMass", "#mu^{+}#mu^{-}K^{#pm}#pi^{#mp} mass [GeV]", "@0*@1 + (1-@0)*@2", RooArgList(tagB0,bMass,bBarMass) )
theBMass     = fulldata.addColumn(theBMassfunc) ;
theBMass.setRange(4.9,5.7);
## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(theBMass,B0Mass) )
deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
deltaPsiPM   = fulldata.addColumn(deltaPMfunc) ;


if args.dimusel == 'keepJpsi':
  cut = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'keepPsiP':
  cut = '(abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( PSIM=PsiPMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'rejectPsi':
  cut = '( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
           (( mumuMass < {JPSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.16 || abs(deltaB0M - deltaPsiPM) < 0.06) ) || \
            ( mumuMass > {PSIM}  && !( abs(deltaB0M - deltaJpsiM) < 0.06 || abs(deltaB0M - deltaPsiPM) < 0.03) ) || \
            ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.06 || abs(deltaB0M - deltaPsiPM) < 0.06 ))))'.format(JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
elif args.dimusel == 'keepPsi':
  cut = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE || abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, PSIM=PsiPMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'nocut':
  cut = 'mumuMass > 0'
else:
  print '\nYou should define which dimuon mass to consider. Please choose between following options: \nkeepPsiP, keepJpsi, rejectPsi, keepPsi'
  sys.exit(0)

print cut


data       = fulldata.reduce(RooArgSet(theBMass,mumuMass,mumuMassE), cut)

mean        = RooRealVar ("mass"         , "mean"          ,  B0Mass_,   3,    7, "GeV")
sigma       = RooRealVar ("#sigma"       , "sigma"         ,  0.028,     0,   10, "GeV")
signalGauss = RooGaussian("signalGauss"  , "signal gauss"  ,  theBMass,  mean,sigma)

sigma2       = RooRealVar ("#sigma2"       , "sigma2"         ,  0.048,     0,   0.07, "GeV")
signalGauss2 = RooGaussian("signalGauss2"  , "signal gauss2"  ,  theBMass,  mean,sigma2)
f1           = RooRealVar ("f1"            , "f1"             ,  0.8  ,     0.,   1.)
gaus         = RooAddPdf  ("gaus"          , "gaus1+gaus2"    , RooArgList(signalGauss,signalGauss2), RooArgList(f1))

pol_c1      = RooRealVar ("p1"           , "coeff x^0 term",    0.5,   -10, 10);
pol_c2      = RooRealVar ("p2"           , "coeff x^1 term",    0.5,   -10, 10);
pol_c3      = RooRealVar ("p3"           , "coeff x^2 term",    0.5,   -10, 10);
slope       = RooRealVar ("slope"        , "slope"         ,    0.5,   -10, 10);
bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   theBMass  );
bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  theBMass, RooArgList(pol_c1,pol_c2));

nsig        = RooRealVar("Yield"         , "signal frac"   ,   4000,     0,   10000);
nbkg        = RooRealVar("nbkg"          , "bkg fraction"  ,   1000,     0,   55000);

fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(signalGauss, bkg_pol), RooArgList(nsig, nbkg))

r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6))
# r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save())#, RooFit.Range(B0Mass_-0.28, B0Mass_+0.28))


frame = theBMass.frame()
data.plotOn(frame, RooFit.Binning(75), RooFit.MarkerSize(.5))
fitFunction.plotOn(frame, );
fitFunction.plotOn(frame, RooFit.Components("bkg_pol")    , RooFit.LineStyle(ROOT.kDashed));
fitFunction.plotOn(frame, RooFit.Components("signalGauss"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));

parList = RooArgSet (nsig,sigma,mean)
# fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+2));

fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
  
frame.GetYaxis().SetTitleOffset(1.35)
frame.getAttText().SetTextSize(0.022) 
frame.getAttText().SetTextFont(42) 
frame.getAttLine().SetLineColor(0) 

c1 = ROOT.TCanvas()
frame.Draw()
c1.SaveAs('save_fit.pdf')

# if args.doubleg==1:
#   resSigma1 = sigma.getVal()
#   resSigma2 = sigma2.getVal()
#   resF      = f1.getVal()
#   
#   totSigma = math.sqrt(resF*(resSigma1**2) + (1-resF)*(resSigma2**2));
#   print 'overall sigma: ', totSigma 


muframe = mumuMass.frame()
data.plotOn(muframe, RooFit.Binning(200), RooFit.MarkerSize(.5))
muframe.Draw()
# c1.SaveAs('save_mumu_2016.pdf')



#pragma link C++ class RooTreeData ;
#pragma link C++ class RooTreeData::PlotOpt ;
#pragma link C++ class RooTruthModel ;
