import argparse

parser = argparse.ArgumentParser(description="")
# parser.add_argument("inputfile" , help = "Path to the input ROOT file")
parser.add_argument("dimusel"   , help = "Define if keep or remove dimuon resonances. You can choose: keepPsiP, keepJpsi, rejectPsi, keepPsi")
parser.add_argument("-d", "--doubleg"  , dest = "doubleg", help = "Define if use 1 or 2 gaussian for the signal model", default = '0')
parser.add_argument("-b", "--use_bdt"  , dest = "use_bdt", help = "Define if use BDT (1) or cut n count (0) selection", default = 1)

args = parser.parse_args()


'''
code to fit the B0 mass distribution:
- unbinned fit
- properly take into account the CP state assignment when plotting the mass variable
- possibility to apply cuts on the dimuon mass [B0&Psi cut in RunI analysis] (e.g. to exclude the Jpsi mass region, or the psi) via the parameter dimusel
'''

import os, sys
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import ROOT
from ROOT import gSystem

gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar
from ROOT import RooGaussian, RooExponential, RooChebychev, RooProdPdf, RooCBShape, TFile, RooPolynomial
import sys
import math	
from uncertainties import ufloat

B0Mass_   = 5.27958
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

B0Mass     = RooRealVar("B0Mass"    , "B0Mass"  , B0Mass_  )
JPsiMass   = RooRealVar("JPsiMass"  , "JPsiMass", JPsiMass_)
PsiPMass   = RooRealVar("PsiPMass"  , "PsiPMass", PsiPMass_)
KStMass    = RooRealVar("KStMass"   , "KStMass" , KStMass_ )

nSigma_psiRej = 3.

q2binning = [
                1,
                2, 
                4.3,
                6,
                8.68,
                10.09,
                12.86,
                14.18,
                16,
#                 19,
]



def fitMC(fulldata, correctTag, ibin):

    print 'now fitting: ', ibin, ' for ', correctTag*'correctTag ', (1-correctTag)*'wrongTag'  
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])

    mean         = RooRealVar ("mass"          , "mean"           ,  B0Mass_,      3,    7, "GeV")
    sigma        = RooRealVar ("#sigma_{1}"    , "sigma"          ,  0.028,        0,   10, "GeV")
    signalGauss  = RooGaussian("signalGauss"   , "signal gauss"   ,  tagged_mass,  mean,sigma)
    
    sigma2       = RooRealVar ("#sigma_{2}"    , "sigma2"         ,  0.048,     0,   0.12, "GeV")
    signalGauss2 = RooGaussian("signalGauss2"  , "signal gauss2"  ,  tagged_mass,  mean,sigma2)
    f1           = RooRealVar ("f1"            , "f1"             ,  0.8  ,     0.,   1.)

    sigma3       = RooRealVar ("#sigma_{3}"    , "sigma3"         ,  0.078,     0,   0.2, "GeV")
    signalGauss3 = RooGaussian("signalGauss3"  , "signal gauss3"  ,  tagged_mass,  mean,sigma3)
    f2           = RooRealVar ("f2"            , "f2"             ,  0.1  ,     0.,   1.)

    sigmaCB      = RooRealVar ("#sigma_{CB}"   , "sigmaCB"        , 0.06  ,     0,   10 )
    alpha        = RooRealVar ("#alpha"        , "alpha"          ,   0.7 ,     0,    5 )
    n            = RooRealVar ("n"	           , "n"	          ,     5 ,     0,   20 )
#     alpha        = RooRealVar ("#alpha"        , "alpha"          ,   2.23,     0,    5 )
#     n            = RooRealVar ("n"	           , "n"	          ,      1,     0,   20 )
    CBShape      = RooCBShape ("CBShape"       , "CBShape"        ,  tagged_mass, mean, sigmaCB, alpha, n)
    f3           = RooRealVar ("f3"            , "f3"             ,  0.3  ,     0.,   1.)

    gaus         = RooAddPdf  ("gaus"          , "gaus1+gaus2"    ,  RooArgList(signalGauss,signalGauss2), RooArgList(f1))
    gausCB       = RooAddPdf  ("gausCB"        , "gaus+CB"        ,  RooArgList(gaus,CBShape), RooArgList(f3))

    pol_c1       = RooRealVar ("p1"            , "coeff x^0 term" ,  -0.5,   -10, 10);
#     bkg_pol      = RooPolynomial("bkg_pol"     , "2nd order pol" ,  tagged_mass);
    bkg_pol     = RooChebychev("bkg_pol"     , "2nd order pol" ,  tagged_mass, RooArgList(pol_c1));
    
    nsig        = RooRealVar("Yield"         , "signal frac"   ,   14000,     0,    1000000)
    nbkg        = RooRealVar("nbkg"          , "bkg fraction"  ,     100,     0,     100000)
    
    
    data        = fulldata.reduce(RooArgSet(thevarsMC), cut)
    if correctTag:
        fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(gaus, bkg_pol), RooArgList(nsig, nbkg))
    else:
        gausCB      = RooAddPdf    ("gausCB"      , "1gaus+CB"      ,  RooArgList(signalGauss,CBShape), RooArgList(f3))
        bkg_pol     = RooPolynomial("bkg_pol"     , "2nd order pol" ,  tagged_mass);
        fitFunction = RooAddPdf    ("fitfunction" , "fit function"  ,  RooArgList(gausCB, bkg_pol), RooArgList(nsig, nbkg))
        
    r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6))
    
    f1_uf    = 0
    s1_uf    = 0
    s2_uf    = 0
    totSigma = 0

    if correctTag and r.status()==0 and r.covQual() == 3:
        f1_uf    = ufloat (f1.getVal(),      f1.getError())
        s1_uf    = ufloat (sigma.getVal(),   sigma.getError())
        s2_uf    = ufloat (sigma2.getVal(),  sigma2.getError())
        totSigma = f1_uf*(s1_uf**2) + (1-f1_uf)*(s2_uf**2)
    if (not correctTag) and r.status()==0 and r.covQual() == 3:
        f3_uf    = ufloat (f3.getVal(),       f3.getError())
        s1_uf    = ufloat (sigma.getVal(),    sigma.getError())
        s2_uf    = ufloat (sigmaCB.getVal(),  sigmaCB.getError())
        totSigma = f3_uf*(s1_uf**2) + (1-f3_uf)*(s2_uf**2)


    print 'sara: status', r.status(), '   covQual ', r.covQual()
    if (r.status()>0 or r.covQual() != 3) and correctTag:
        print 'RT: going to retry the fit with 3 gaussians for bin %s'%ibin
        sum3gaus    = RooAddPdf ("gaus3"       , "gaus+gaus3"    ,  RooArgList(signalGauss, signalGauss2, signalGauss3), RooArgList(f1,f2))
        fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(sum3gaus, bkg_pol), RooArgList(nsig, nbkg))
        r           = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6))
        if r.status()==0 and r.covQual() == 3:
            f1_uf    = ufloat (f1.getVal(),      f1.getError())
            f2_uf    = ufloat (f2.getVal(),      f2.getError())
            s1_uf    = ufloat (sigma.getVal(),   sigma.getError())
            s2_uf    = ufloat (sigma2.getVal(),  sigma2.getError())
            s3_uf    = ufloat (sigma3.getVal(),  sigma3.getError())
            totSigma = f1_uf*(s1_uf**2) + f2_uf*s2_uf**2  + (1-f1_uf-s2_uf)*(s3_uf**2)

  
    if (r.status()>0 or r.covQual() != 3) and not correctTag:
        print 'WT: going to retry the fit with 2 gaussian for bin %s'%ibin
        gausCB      = RooAddPdf ("gausCB"      , "gaus+CB"       ,  RooArgList(gaus,CBShape)   , RooArgList(f3))
        fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(gausCB, bkg_pol), RooArgList(nsig, nbkg))
        r           = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6))
    
    frame = tagged_mass.frame()
    data.plotOn(frame, RooFit.Binning(60), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, RooFit.Components("bkg_pol")      , RooFit.LineStyle(ROOT.kDashed));
    fitFunction.plotOn(frame, RooFit.Components("signalGauss")  , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));
    fitFunction.plotOn(frame, RooFit.Components("signalGauss2") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+2));
    color = ROOT.kGreen+2
    fitFunction.plotOn(frame, RooFit.Components("CBShape")      , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+3));

    if correctTag:
        color = ROOT.kMagenta+2

#     fitFunction.plotOn(frame, RooFit.Components("gaus") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(color), RooFit.FillColor(color), RooFit.DrawOption("F"), RooFit.FillStyle(3004))
    fitFunction.plotOn(frame, RooFit.Components("gaus")   , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(color));
    fitFunction.plotOn(frame, RooFit.Components("gausCB") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kAzure+3));

    fitFunction.plotOn(frame )
    
    print 'now paramON: '
#     parList = RooArgSet (nsig,sigma,sigma2, sigmaCB, mean)
    fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88))
#     fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
    frame.Draw()
    
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.getAttText().SetTextSize(0.022) 
    frame.getAttText().SetTextFont(42) 
    frame.getAttLine().SetLineColor(0) 

# #     dict_s_v1[ibdt]  = [nsig.getVal() / mclumi*datalumi , nsig.getError()/ mclumi*datalumi]
    f1_uf    = ufloat (f1.getVal(),      f1.getError())
    s1_uf    = ufloat (sigma.getVal(),   sigma.getError())
    s2_uf    = ufloat (sigma2.getVal(),  sigma2.getError())
    totSigma = f1_uf*(s1_uf**2) + (1-f1_uf)*(s2_uf**2)
# 
    if correctTag:
        dict_s_rt[ibin]   = [nsig.getVal(), nsig.getError()]
        sigma_rt_mc[ibin] = [ math.sqrt(totSigma.n), totSigma.s/2/math.sqrt(totSigma.n)]
        frame.SetTitle('correctly tagged events')
        c1.SaveAs('fit_results_mass/save_fit_mc_%s_2018_RT.pdf'%ibin)
    else:
        dict_s_wt[ibin]     = [nsig.getVal(), nsig.getError()]
        sigma_wt_mc[ibin] = [ math.sqrt(totSigma.n), totSigma.s/2/math.sqrt(totSigma.n)]
        frame.SetTitle('wrongly tagged events')
        c1.SaveAs('fit_results_mass/save_fit_mc_%s_2018_WT.pdf'%ibin)

    out_f.cd()
    r.Write('results_%s_%s'%(correctTag*'RT' + (1-correctTag)*'WT', ibin))

   
def fitData(fulldata, ibin):

#     import pdb; pdb.set_trace()
    n_rt = ufloat(dict_s_rt[ibin][0], dict_s_rt[ibin][1])
    n_wt = ufloat(dict_s_wt[ibin][0], dict_s_wt[ibin][1])
    fraction = n_rt/(n_rt+n_wt)
    print 'mistag fraction on MC for bin ', ibin , ' : ' , fraction.n , fraction.s 
    
    mean          = RooRealVar ("mass"          , "mean"            ,  B0Mass_,   4,   6, "GeV")

    sigmart       = RooRealVar ("#sigma_{rt}"   , "sigmart"         , sigma_rt_mc[ibin][0], 0, 1.)
    c_sigma_rt    = RooGaussian("c_sigma_rt"    , "c_sigma_rt"      , sigmart, ROOT.RooFit.RooConst(sigma_rt_mc[ibin][0]), ROOT.RooFit.RooConst(sigma_rt_mc[ibin][1])) ;
    rtGauss       = RooGaussian("rtGauss"       , "signal gauss"    , tagged_mass,  mean, sigmart)
    c_rtGauss     = RooProdPdf ("c_rtGauss"     , "constr rt gauss" , rtGauss,      c_sigma_rt )     

    sigmawt       = RooRealVar  ("#sigma_{wt}"  , "sigmawt"         , sigma_wt_mc[ibin][0], 0, 1)
    c_sigma_wt    = RooGaussian ("c_sigma_wt"   , "c_sigma_wt"      , sigmawt,  ROOT.RooFit.RooConst(sigma_wt_mc[ibin][0]), ROOT.RooFit.RooConst(sigma_wt_mc[ibin][1])) ;
    wtGauss       = RooGaussian ("wtGauss"      , "signalGausswt"   , tagged_mass,  mean, sigmawt )
    c_wtGauss     = RooProdPdf  ("c_wtGauss"    , "constr wt gauss" , wtGauss, c_sigma_wt     )     

    f1            = RooRealVar  ("f1"           , "f1"              , fraction.n, 0, 1)
    c_f1          = RooGaussian ("c_f1"         , "c_f1"            , f1,  ROOT.RooFit.RooConst(fraction.n ), ROOT.RooFit.RooConst(fraction.s)) ;
    gaus          = RooAddPdf   ("gaus"         , "gaus1+gaus2"     , RooArgList(c_rtGauss,c_wtGauss), RooArgList(f1))
    c_gaus        = RooProdPdf  ("c_gaus"       , "c_gaus"          , gaus, c_f1)
    
    slope         = RooRealVar ("slope"         , "slope"           ,    0.5,   -10, 10);
    bkg_exp       = RooExponential("bkg_exp"    , "exponential"     ,  slope,   tagged_mass  );
    pol_c1        = RooRealVar  ("p1"           , "coeff x^0 term"  ,    0.5,   -10, 10);
    pol_c2        = RooRealVar  ("p2"           , "coeff x^1 term"  ,    0.5,   -10, 10);
    bkg_pol       = RooChebychev("bkg_pol"      , "2nd order pol"   ,  tagged_mass, RooArgList(pol_c1,pol_c2));
   
    nsig         = RooRealVar("Yield"         , "signal frac"    ,    4000,     0,   1000000);
    nbkg         = RooRealVar("nbkg"          , "bkg fraction"   ,    1000,     0,   550000);
    
    cut = cut_base + '&& (mumuMass*mumuMass > %s && mumuMass*mumuMass < %s)'%(q2binning[ibin], q2binning[ibin+1])
    
    data       = fulldata.reduce(RooArgSet(tagged_mass,mumuMass,mumuMassE), cut)

    fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(c_gaus, bkg_exp), RooArgList(nsig, nbkg))

    r = fitFunction.fitTo(data, 
                          RooFit.Extended(True), 
                          RooFit.Save(), 
                          RooFit.Range(5.,5.6), 
                          RooFit.Verbose(True),
#                           ROOT.RooFit.Constrain(RooArgSet(sigmart,f1))
                          ROOT.RooFit.Constrain(RooArgSet(sigmart,sigmawt,f1))
                         )
#     RooFitResult* r2 = modelc.fitTo(*d,Constrain(f),Save()) ;
    
    frame = tagged_mass.frame()
    data.plotOn(frame, RooFit.Binning(60), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame);
#     fitFunction.plotOn(frame, RooFit.Components("signalGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1))
    fitFunction.plotOn(frame, RooFit.Components("rtGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1));
    fitFunction.plotOn(frame, RooFit.Components("wtGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kGreen+1))
    fitFunction.plotOn(frame, RooFit.Components("gaus"   ) , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kBlue), RooFit.DrawOption("F"), RooFit.FillStyle(3004))
    fitFunction.plotOn(frame, RooFit.Components("gaus"   ) , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kBlue));
    fitFunction.plotOn(frame, RooFit.Components("bkg_exp") , RooFit.LineStyle(ROOT.kDashed));
#     fitFunction.plotOn(frame, RooFit.Components("signalGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1));
#     fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1));
    
    parList = RooArgSet (nsig,sigmart,sigmawt, mean)
    # fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+2));
    
    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
    frame.Draw()
      
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.getAttText().SetTextSize(0.022) 
    frame.getAttText().SetTextFont(42) 
    frame.getAttLine().SetLineColor(0) 
    frame.SetTitle('')

    txt = ROOT.TLatex(.11,.91,"CMS") ;
    txt. SetNDC() ;
    txt. SetTextSize(0.045) ;
    frame.addObject(txt) ;
        
    txt2 = ROOT.TLatex(.75,.91,"61.1 fb^{-1}, 13 TeV") ;
    txt2 . SetNDC() ;
    txt2 . SetTextSize(0.03) ;
    txt2 . SetTextFont(42) ;
    frame. addObject(txt2) ;

#     txtq = ROOT.TLatex(.7, .41, '%s' %ibin) ;
    txtq = ROOT.TLatex(.6,.5, "%s GeV^{2} < q^{2} < %s GeV^{2}" %(q2binning[ibin], q2binning[ibin+1])) ;
    txtq . SetNDC() ;
    txtq . SetTextSize(0.033) ;
    txtq . SetTextFont(42) ;
    frame. addObject(txtq) ;
       
    frame.Draw()
    c1.SaveAs('fit_results_mass/save_fit_data_%s_2018_LMNR_c1.pdf'%ibin)
 
    resSigma1 = ufloat (sigmart.getVal(),  sigmart.getError())
    resSigma2 = ufloat (sigmawt.getVal(), sigmawt.getError())
    resF      = ufloat (f1.getVal(),     f1.getError())

    totSigma = resF*(resSigma1**2) + (1-resF)*(resSigma2**2)

    totSigma_v = math.sqrt(totSigma.n)
    totSigma_e = totSigma.s/2/math.sqrt(totSigma.n)
 
    yields[ibin] =  [nsig.getVal(), nsig.getError()]
    sigmas[ibin] =  [totSigma_v   , totSigma_e     ]
    f1s[ibin]    =  [f1.getVal()  , f1.getError()  ]
    
    sigma_rt_data[ibin] = [sigmart.getVal()  , sigmart.getError()]
    sigma_wt_data[ibin] = [sigmawt.getVal(), sigmawt.getError()]




tData = ROOT.TChain('ntuple')
tData.Add('final_ntuples/2018Data_All_finalSelection.root')
## test sample
# tData.Add('final_ntuples/2018Data_All_finalSelection_100000.root')

tMC = ROOT.TChain('ntuple')
tMC.Add('final_ntuples/2018MC_LMNR.root')
## test sample
# tMC.Add('final_ntuples/2018MC_LMNR_5000.root')

tagged_mass     = RooRealVar("tagged_mass"    , "#mu^{+}#mu^{-}K* mass", 4.9, 5.6, "GeV")
mumuMass        = RooRealVar("mumuMass" , "mumuMass" , 0, 6);
mumuMassE       = RooRealVar("mumuMassE", "mumuMassE", 0, 10000);
tagB0           = RooRealVar("tagB0"    , "tagB0"    , 0, 2);

thevars = RooArgSet()
thevars.add(tagged_mass)
thevars.add(mumuMass)
thevars.add(mumuMassE)
thevars.add(tagB0)

fulldata   = RooDataSet('fulldata', 'fulldataset', tData,  RooArgSet(thevars))


## add to the input tree the combination of the variables for the B0 arb. mass
# theBMassfunc = RooFormulaVar("theBMass", "#mu^{+}#mu^{-}K^{#pm}#pi^{#mp} mass", "@0*@1 + (1-@0)*@2", RooArgList(tagB0,tagged_mass) )
# theBMass     = fulldata.addColumn(theBMassfunc) ;
# theBMass.setRange(5.0,5.6);
# theBMass.setUnit("GeV")
## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(tagged_mass,B0Mass) )
deltaB0M     = fulldata.addColumn(deltaB0Mfunc) ;
deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
deltaJpsiM   = fulldata.addColumn(deltaJMfunc) ;
deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
deltaPsiPM   = fulldata.addColumn(deltaPMfunc) ;

genSignal   = RooRealVar("genSignal"    , "genSignal"    , 0, 10);
thevarsMC   = thevars; thevarsMC.add(genSignal)
fullmc      = RooDataSet('fullmc', 'fullmc', tMC,  RooArgSet(thevarsMC))
# tagged_mass    = fullmc.addColumn(theBMassfunc) 
deltaB0M    = fullmc.addColumn(deltaB0Mfunc) 
deltaJpsiM  = fullmc.addColumn(deltaJMfunc)  
deltaPsiPM  = fullmc.addColumn(deltaPMfunc)  

thevars.add(deltaB0M)
thevars.add(deltaJpsiM)
thevars.add(deltaPsiPM)

thevarsMC.add(deltaB0M)
thevarsMC.add(deltaJpsiM)
thevarsMC.add(deltaPsiPM)


rt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '(tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2)')
wt_mc       = fullmc.reduce(RooArgSet(thevarsMC), '(tagB0==0 && genSignal==1) || (tagB0==1 && genSignal==2)')

# import pdb; pdb.set_trace()

if args.dimusel == 'keepJpsi':
  cut_base = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'keepPsiP':
  cut_base = '(abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( PSIM=PsiPMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'rejectPsi':
  cut_base = '( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
           (( mumuMass < {JPSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.18 || abs(deltaB0M - deltaPsiPM) < 0.0) ) || \
            ( mumuMass > {PSIM}  && !( abs(deltaB0M - deltaJpsiM) < 0.0  || abs(deltaB0M - deltaPsiPM) < 0.09) ) || \
            ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.08 || abs(deltaB0M - deltaPsiPM) < 0.08 ))))'.format(JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
elif args.dimusel == 'keepPsi':
  cut_base = '(abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE || abs(mumuMass - {PSIM}) < {CUT}*mumuMassE)'.format( JPSIM=JPsiMass_, PSIM=PsiPMass_, CUT=nSigma_psiRej)
elif args.dimusel == 'nocut':
  cut_base = 'mumuMass > 0'
else:
  print '\nYou should define which dimuon mass to consider. Please choose between following options: \nkeepPsiP, keepJpsi, rejectPsi, keepPsi'
  sys.exit(0)

# cut = cut + '&& (mumuMass*mumuMass > 16 && mumuMass*mumuMass < 19)'



c1 = ROOT.TCanvas() 

yields = {}
sigmas = {}
f1s    = {}
sigma_rt_data = {}
sigma_wt_data = {}

dict_s_rt  = {}
dict_s_wt  = {}
sigma_rt_mc = {}
sigma_wt_mc = {}

out_f = TFile ("results_fits.root","RECREATE") 

for ibin in range(len(q2binning)-1):


    print 'dimuon selection: ', args.dimusel
    if args.dimusel == 'rejectPsi' and \
       (q2binning[ibin] == 8.68 or q2binning[ibin] == 12.86): 
           continue
           
    fitMC(rt_mc, True, ibin)
    fitMC(wt_mc, False, ibin)

    fitData(fulldata, ibin)



 
 
for k,v in yields.items():
    print 'bin ', k, '\t B0 yield: ', v[0], '+/-', v[1]

for k,v in sigmas.items():
    print 'bin ', k, '\t B0 sigma: ', v[0], '+/-', v[1]

for k,v in sigma_rt_mc.items():
    print 'bin ', k, '\t B0 sigma MC rt : ', v[0], '+/-', v[1]
for k,v in sigma_rt_data.items():
    print 'bin ', k, '\t B0 sigma data rt : ', v[0], '+/-', v[1]

for k,v in sigma_wt_mc.items():
    print 'bin ', k, '\t B0 sigma MC wt : ', v[0], '+/-', v[1]
for k,v in sigma_wt_data.items():
    print 'bin ', k, '\t B0 sigma data wt : ', v[0], '+/-', v[1]

for i in range(len(dict_s_rt)):
    
    n_rt = ufloat(dict_s_rt[ibin][0], dict_s_rt[ibin][1])
    n_wt = ufloat(dict_s_wt[ibin][0], dict_s_wt[ibin][1])
    fraction = n_rt/(n_rt+n_wt)
    print 'bin ', ibin, '\t mistag fraction MC '  , ibin , ' : ' , fraction.n , '+/-', fraction.s 
    print 'bin ', ibin, '\t mistag fraction data ', ibin , ' : ' , f1s[i][0], '+/-', f1s[i][1]



out_f.Close() 


# if args.doubleg==1:
#   resSigma1 = sigma.getVal()
#   resSigma2 = sigma2.getVal()
#   resF      = f1.getVal()
#   
#   totSigma = math.sqrt(resF*(resSigma1**2) + (1-resF)*(resSigma2**2));
#   print 'overall sigma: ', totSigma 


# muframe = mumuMass.frame()
# data.plotOn(muframe, RooFit.Binning(200), RooFit.MarkerSize(.5))
# muframe.Draw()
# c1.SaveAs('save_mumu_2016.pdf')



#pragma link C++ class RooTreeData ;
#pragma link C++ class RooTreeData::PlotOpt ;
#pragma link C++ class RooTruthModel ;
