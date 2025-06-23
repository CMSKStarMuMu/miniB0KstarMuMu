import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("year" , help = "year")
parser.add_argument("-c", "--channel" , dest = "channel",  help = "Define if Jpsi or LMNR or PsiPrime", default='Jpsi')
parser.add_argument("-e", "--era"     , dest = "era"    ,  help = "Define BH, GH, BF"                 , default='BH')
parser.add_argument("-l", "--l1"      , dest = "l1seed"     , help = "L1 seed: l1_11_4, l1_12_5, l1_10_0, l1_00, l1_00_OS", default = 'all')

args = parser.parse_args()

import ROOT
import sys
import math, pdb
from array import array


B0Mass_   = 5.27963
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

b0_width   = 0.03601 
kaonMass_  = 0.493677

mc_sigma = 0.040
mc_mass  = 5.27783 

nSigma_psiRej = 3.

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)


class PlotContainer(object):
    def __init__(self, var, sel, nanosel, xtitle, norm, pltopt, nbins=0, xmin=0, xmax=0, mybins=0, label=None, logx=False, logy=False, fill=False):
        self.var         = var 
        self.sel         = sel 
        self.nanosel       = nanosel 
        self.xtitle      = xtitle 
        self.norm        = norm 
        self.nbins       = nbins 
        self.xmin        = xmin 
        self.xmax        = xmax
        self.pltopt      = pltopt
        self.logy        = logy
        self.logx        = logx
        self.fill        = fill
        self.mybins      = mybins 
        if label:
            self.label = label
        else:
            self.label = self.var

def setHistoProperties(h, color, xtitle, ytitle, xoff, yoff, markerSt, markerSize ):
    h.SetLineColor(color )
    h.SetFillColor(color )
    h.SetMarkerColor(color)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.GetXaxis().SetTitleOffset(xoff)
    h.GetYaxis().SetTitleOffset(yoff)
    h.SetMarkerStyle(markerSt)
    h.SetMarkerSize(markerSize)
    h.SetMarkerColor(color)
    h.GetXaxis().SetLabelColor(ROOT.kWhite)
    h.GetXaxis().SetLabelSize(0.)
    h.SetNdivisions(0)

def create_histo(name, nbins, xmin, xmax, mybins):
    h_data = ROOT.TH1F(name+'_%d'%i    , '', nbins, xmin, xmax)
    if mybins == 0:
      hr_data  = ROOT.TH1F(name+'_ratio_%d'%i   , '', nbins, xmin, xmax)
    else:
#       h_data = ROOT.TH1F(name+'_%d'%i    , '', len(mybins)-1, array('d',mybins))
      hr_data  = ROOT.TH1F(name+'_ratio_%d'%i   , '', len(mybins)-1, array('d',mybins))
            
    return h_data, hr_data


tnano  = ROOT.TChain('ntuple')
tntuple  = ROOT.TChain('ntuple') 

tnano .Add('/gwpool/users/fiorendi/p5prime/run3/CMSSW_14_2_2/src/miniB0KstarMuMu/miniKstarMuMu/test/flat_ntuples/ntuple_flat.root')
tntuple .Add('ntuple_flat_ourNtuple.root')

print ('nano', tnano.GetName()  )
print ('ntuple', tntuple.GetName()  )

pt_mu1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(10)]
pt_mu2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(10)] + [i*2+20 for i in range(1)] + [i*4+22 for i in range(3)]
pt_tk1_bins = [i*0.1 for i in range(50)]
pt_tk2_bins = [i*0.1 for i in range(50)]
# pt_tk1_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] + [i*2.5+15 for i in range(3)]
# pt_tk2_bins = [i*0.5 for i in range(20)] + [i+10 for i in range(5) ] 

L_bins     = [i*0.01 for i in range(30)] + [i*0.05+0.3 for i in range(2)] #+ [i*0.1+1. for i in range(6)]
LE_bins    = [0] + [0.002+i*0.0004 for i in range(5)] + [i*0.0002+0.004 for i in range(20)] + [i*0.0004+0.008 for i in range(5)] + [i*0.001+0.01 for i in range(4)] + [i*0.002+0.014 for i in range(4)]
LS_bins    = [i*2 for i in range(20)]    #+ [i*4+100 for i in range(5)] #  + [i*10+180. for i in range(3)]
DCA_bins   = [i*0.05 for i in range(50)] + [i*0.1+2.5 for i in range(10)] + [i*0.5+3.5 for i in range(10)]

cos_bins   = [0.89, 0.9, 0.994, 0.997, 0.998, 0.999, 0.9992, 0.9994, 0.9996, 0.9997, 0.9998, 0.9999, 0.99995, 1.]

selNANO   = 'bCosAlphaBS>0.9 && mumPt > 3 && mupPt > 3 && bVtxCL > 0.001 && kstTrkmPt > 0.8 && kstTrkpPt > 0.8 && kstTrkpDCABSSign > 0.8 && kstTrkmDCABSSign > 0.8'
selNtuple = 'bCosAlphaBS>0.9 && mumPt > 3 && mupPt > 3 && bVtxCL > 0.001 && kstTrkmPt > 0.8 && kstTrkpPt > 0.8 && kstTrkpDCABSSign > 0.8 && kstTrkmDCABSSign > 0.8'
    
selNANO += ' && abs(mumEta) < 2.4  && abs(mupEta) < 2.4 && mumSoft > 0 && mupSoft > 0'
selNtuple += ' && abs(mumEta) < 2.4  && abs(mupEta) < 2.4'

# selNANO += '&& eventN == 15496480'
# selNtuple += '&& eventN == 15496480'
# import pdb
# pdb.set_trace()   
# selNtuple = '((tagB0==1 && (bMass    > {M}-2.5*{S} && bMass    < {M}+2.5*{S}) ) || \
#             (tagB0==0 && (bBarMass > {M}-2.5*{S} && bBarMass < {M}+2.5*{S}) ))\
#           '.format(M=mc_mass,S=mc_sigma)
#           (kstMass*tagB0 + kstBarMass*(1-tagB0) > 0.806 && kstMass*tagB0 + kstBarMass*(1-tagB0) < 0.986 )  '.format(M=mc_mass,S=mc_sigma)


thebmass     = 'bMass*tagB0 + bBarMass*(1-tagB0)'
thekstmass   = 'kstMass*tagB0 + kstBarMass*(1-tagB0)'


toplot = [
    PlotContainer(var = 'bPt'                   , sel = selNtuple, nanosel = selNANO, xtitle = 'p_{T}(B^{0}) [GeV]'           , norm = True, nbins = 150 , xmin =  0            , xmax =  80  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bEta'                  , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(B^{0})'                  , norm = True, nbins = 50 , xmin =  -2.5         , xmax =  2.5 , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bVtxCL'                , sel = selNtuple, nanosel = selNANO, xtitle = 'B^{0} vtx CL'                 , norm = True, nbins = 100  , xmin =  0            , xmax =  1   , fill = True, pltopt = 'HIST', logy = True, logx = True),
    PlotContainer(var = 'kstVtxCL'              , sel = selNtuple, nanosel = selNANO, xtitle = 'K* vtx CL'                    , norm = True, nbins = 100  , xmin =  0            , xmax =  1   , fill = True, pltopt = 'HIST', logy = True, logx = True),
    PlotContainer(var = 'mumuVtxCL'             , sel = selNtuple, nanosel = selNANO, xtitle = 'dimuon vtx CL'                , norm = True, nbins = 100  , xmin =  0            , xmax =  1   , fill = True, pltopt = 'HIST', logy = True, logx = True),
    PlotContainer(var = 'mumPt'                 , sel = selNtuple, nanosel = selNANO, xtitle = 'p_{T}(#mu_{-}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu1_bins[0], xmax =  pt_mu1_bins[-1] , mybins = pt_mu1_bins, fill = True, pltopt = 'HIST', label = 'leadingMuPt' ),
    PlotContainer(var = 'mupPt'                 , sel = selNtuple, nanosel = selNANO, xtitle = 'p_{T}(#mu_{+}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu2_bins[0], xmax =  pt_mu2_bins[-1] , mybins = pt_mu2_bins, fill = True, pltopt = 'HIST', label = 'trailingMuPt'),
    PlotContainer(var = 'kstTrkmPt'             , sel = selNtuple, nanosel = selNANO, xtitle = 'trk- p_{T} [GeV]'             , norm = True, nbins = 80  , xmin = pt_tk1_bins[0], xmax =  pt_tk1_bins[-1] , mybins = pt_tk1_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpPt'             , sel = selNtuple, nanosel = selNANO, xtitle = 'trk+ p_{T} [GeV]'             , norm = True, nbins = 80  , xmin = pt_tk2_bins[0], xmax =  pt_tk2_bins[-1] , mybins = pt_tk2_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumEta'                , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(#mu_{-})'                , norm = True, nbins = 50 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingMuEta' ),
    PlotContainer(var = 'mupEta'                , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(#mu_{+})'                , norm = True, nbins = 50 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingMuEta'),
    PlotContainer(var = 'kstTrkmEta'            , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(tk_{-})'                 , norm = True, nbins = 50 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingTkEta' ),
    PlotContainer(var = 'kstTrkpEta'            , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(tk_{+})'                 , norm = True, nbins = 50 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingTkEta'),
#     PlotContainer(var = 'mu1Pt'                 , sel = selNtuple, nanosel = selNANO, xtitle = 'p_{T}(#mu_{1}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu1_bins[0], xmax =  pt_mu1_bins[-1] , mybins = pt_mu1_bins, fill = True, pltopt = 'HIST', label = 'leadingMuPt' ),
#     PlotContainer(var = 'mu2Pt'                 , sel = selNtuple, nanosel = selNANO, xtitle = 'p_{T}(#mu_{2}) [GeV]'         , norm = True, nbins = 160 , xmin = pt_mu2_bins[0], xmax =  pt_mu2_bins[-1] , mybins = pt_mu2_bins, fill = True, pltopt = 'HIST', label = 'trailingMuPt'),
#     PlotContainer(var = 'kstTrk1Pt'             , sel = selNtuple, nanosel = selNANO, xtitle = 'leading trk p_{T} [GeV]'      , norm = True, nbins = 80  , xmin = pt_tk1_bins[0], xmax =  pt_tk1_bins[-1] , mybins = pt_tk1_bins,   fill = True, pltopt = 'HIST'),
#     PlotContainer(var = 'kstTrk2Pt'             , sel = selNtuple, nanosel = selNANO, xtitle = 'trailing trk p_{T} [GeV]'     , norm = True, nbins = 80  , xmin = pt_tk2_bins[0], xmax =  pt_tk2_bins[-1] , mybins = pt_tk2_bins,   fill = True, pltopt = 'HIST'),
#     PlotContainer(var = 'mu1Eta'                , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(#mu_{1})'                , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingMuEta' ),
#     PlotContainer(var = 'mu2Eta'                , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(#mu_{2})'                , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingMuEta'),
#     PlotContainer(var = 'kstTrk1Eta'            , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(tk_{1})'                 , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'leadingTkEta' ),
#     PlotContainer(var = 'kstTrk2Eta'            , sel = selNtuple, nanosel = selNANO, xtitle = '#eta(tk_{2})'                 , norm = True, nbins = 100 , xmin =  -2.5        , xmax =  2.5 , fill = True, pltopt = 'HIST', label = 'trailingTkEta'),
    PlotContainer(var = 'bCosAlphaBS'           , sel = selNtuple, nanosel = selNANO, xtitle = 'cos#theta_{BS}'               , norm = True, nbins = 100 , xmin = cos_bins[0]  , xmax = cos_bins[-1], mybins = cos_bins, fill = True, pltopt = 'HIST', logy = True, label = 'bCosThetaBS'),
    PlotContainer(var = 'bLBS'                  , sel = selNtuple, nanosel = selNANO, xtitle = 'L_{BS} [cm]'                  , norm = True, nbins = 150 , xmin = L_bins[0]    , xmax = L_bins[-1]  , mybins = L_bins,   fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bLBSE'                 , sel = selNtuple, nanosel = selNANO, xtitle = 'L_{BS} uncertainty [cm]'      , norm = True, nbins = 100 , xmin = LE_bins[0]   , xmax = LE_bins[-1] , mybins = LE_bins,  fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bLBS/bLBSE'            , sel = selNtuple, nanosel = selNANO, xtitle = 'L_{BS}/#sigma_{L}'            , norm = True, nbins = 100 , xmin = LS_bins[0]   , xmax = LS_bins[-1] , mybins = LS_bins,  fill = True, pltopt = 'HIST', label = 'LSigmaBS'),
    PlotContainer(var = 'bDCABS/bDCABSE'        , sel = selNtuple, nanosel = selNANO, xtitle = 'b DCA from BS significance'   , norm = True, nbins = 100 , xmin = DCA_bins[0]  , xmax = DCA_bins[-1], mybins = DCA_bins, fill = True, pltopt = 'HIST', label = 'bDCASign'),
    PlotContainer(var = 'mumuMass'              , sel = selNtuple, nanosel = selNANO, xtitle = 'mumuMass'                     , norm = True, nbins = 100 , xmin = 0  , xmax = 5, fill = True, pltopt = 'HIST', label = 'mumuMass', logy = True),
    PlotContainer(var = thekstmass              , sel = selNtuple, nanosel = selNANO, xtitle = 'kstar mass [GeV]'             , norm = True, nbins = 120  , xmin =  0.4         , xmax =  1.6 , fill = True, pltopt = 'HIST', label = 'kstMass'),
    PlotContainer(var = thebmass                , sel = selNtuple, nanosel = selNANO, xtitle = 'b mass [GeV]'                 , norm = True, nbins = 100  , xmin =  4.         , xmax =  6 , fill = True, pltopt = 'HIST', label = 'bMass'),
#     PlotContainer(var = 'sum_isopt_04'          , sel = selNtuple, nanosel = selNANO, xtitle = 'B^{0} isolation'              , norm = True, nbins = 50  , xmin =  0           , xmax =  5   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpDCABSSign'      , sel = selNtuple, nanosel = selNANO, xtitle = 'trk+ DCA(BS) sign [cm]'        , norm = True, nbins = 80  , xmin =  -0.4         , xmax =  10            ,                    fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmDCABSSign'      , sel = selNtuple, nanosel = selNANO, xtitle = 'trk- DCA(BS) sign [cm]'        , norm = True, nbins = 80  , xmin =  -0.4         , xmax =  10            ,                    fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpDCABS'          , sel = selNtuple, nanosel = selNANO, xtitle = 'trk+ DCA(BS) [cm]'            , norm = True, nbins = 100  , xmin =  -1         , xmax =  1            ,                    fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmDCABS'          , sel = selNtuple, nanosel = selNANO, xtitle = 'trk- DCA(BS) [cm]'            , norm = True, nbins = 100  , xmin =  -1         , xmax =  1            ,                    fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kkMass'                , sel = selNtuple, nanosel = selNANO, xtitle = 'KK mass [GeV]'               , norm = True, nbins = 90  , xmin =  0.6         , xmax =  2.4 , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'dR_mum_trkm'           , sel = selNtuple, nanosel = selNANO, xtitle = 'dR_mum_trkm'               , norm = True, nbins = 200  , xmin =  0.         , xmax =  2.5 , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'cos_theta_l'           , sel = selNtuple, nanosel = selNANO, xtitle = 'cos_theta_l'                , norm = True, nbins = 100  , xmin =  -1            , xmax =  1   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'runN'                  , sel = selNtuple, nanosel = selNANO, xtitle = 'runN'                         , norm = True, nbins = 50  , xmin =  392150      , xmax = 392200         ,                    fill = True, pltopt = 'HIST'),

]



ROOT.TH1.SetDefaultSumw2()
colorlist = [ROOT.kBlack,  ROOT.kMagenta+2, ROOT.kOrange-3, ROOT.kGreen+3,ROOT.kCyan-2, ROOT.kBlue-3, ROOT.kViolet+5, ROOT.kSpring-5]


for i, iplot in enumerate(toplot):
    
    nbins  = iplot.nbins
    xmin   = iplot.xmin  
    xmax   = iplot.xmax  
    sel    = iplot.sel
    nanosel  = iplot.nanosel
    var    = iplot.var
    
    h_nano, hr_nano = create_histo('hnano', nbins, xmin, xmax, iplot.mybins)
    h_ntuple, hr_ntuple = create_histo('hntuple', nbins, xmin, xmax, iplot.mybins)

    xtitle = iplot.xtitle
    ytitle = 'a.u.' if iplot.norm else 'events'
    
    p_th1  = ROOT. TH1F("rp_th1","rp_th1",nbins,xmin,xmax)
    setHistoProperties(h_nano    , colorlist[0], xtitle, ytitle, 1.2, 1.5, 20, 0.6 )
    setHistoProperties(h_ntuple  , colorlist[1], xtitle, ytitle, 1.2, 1.5, 22, 0.6 )
    
    tnano.Draw('%s >> %s' %(var, h_nano.GetName())      , '(%s)'   %(nanosel)   , iplot.pltopt   )
    tntuple.Draw('%s >> %s' %(var, h_ntuple.GetName())    , '(%s)'   %(sel)   , iplot.pltopt   )

    tnano.Draw('%s >> %s' %(var, hr_nano.GetName())     , '(%s)'   %(nanosel)              , iplot.pltopt   )
    tntuple.Draw('%s >> %s' %(var, hr_ntuple.GetName())   , '(%s)'   %(sel)              , iplot.pltopt   )

    ymax = max(h_nano.GetMaximum(), h_ntuple.GetMaximum())
    p_th1 .SetTitle("");
    p_th1.GetYaxis().SetRangeUser(0.01,1.2*ymax)
    p_th1.SetMaximum(1.2*ymax)
    p_th1.GetXaxis().SetLabelColor(ROOT.kWhite)

       
    c1 = ROOT.TCanvas('c1', 'c1', 700, 700)
    upperPad  = ROOT.TPad('upperPad' , 'upperPad' , 0., 0.35 , 1.,  1.    )  
    lowerPad  = ROOT.TPad('lowerPad' , 'lowerPad' , 0., 0.0  , 1.,  0.345 )  
    upperPad.Draw()
    lowerPad .Draw()
    upperPad.SetBottomMargin(0.012)
    lowerPad.SetTopMargin(0)
    lowerPad.SetBottomMargin(0.2)
   
    upperPad.cd()
    p_th1.Draw()

    h_nano .Draw('e same')
    h_ntuple .Draw('e same')


    l1 = ROOT.TLegend(0.55,0.73,0.88,0.88)
    l1.AddEntry(h_nano, 'nano' , 'p')
    l1.AddEntry(h_ntuple, 'ntuple' , 'p')
    l1.SetBorderSize(0)
    l1.SetTextSize(0.036)
    l1.Draw('same')
    
    if iplot.logy:
      upperPad.SetLogy()
    if iplot.logx:
      upperPad.SetLogx()

    c1.Update()
    c1.Modified()


    r_nano    = ROOT.TGraphAsymmErrors(hr_nano.GetNbinsX())
    r_ntuple  = ROOT.TGraphAsymmErrors(hr_ntuple.GetNbinsX())

    hr_data_mc_b = hr_nano.Clone()
    hr_data_mc_b.Divide(hr_nano, hr_ntuple, 1, 1, 'B')
    r_b = ROOT.TGraphAsymmErrors(hr_data_mc_b)

    r_b .SetMarkerStyle(8)
    r_b .SetMarkerColor(colorlist[1])
    r_b .SetLineColor  (colorlist[1])

    rp_th1  = ROOT. TH1F("rp_th1","rp_th1",nbins,xmin,xmax)
    rp_th1 .SetTitle("");
    rp_th1 .SetDirectory(0);
    rp_th1 .SetStats(0);
    rp_th1 .GetXaxis().SetTitle(xtitle);
    rp_th1 .GetYaxis().SetTitle('nano/ntuple');
    rp_th1 .GetXaxis().SetTitleSize(0.08);
    rp_th1 .GetYaxis().SetTitleSize(0.07);
    rp_th1 .GetYaxis().SetTitleOffset(0.7);
    rp_th1 .GetYaxis().SetNdivisions(4+800)
    rp_th1 .GetXaxis().SetTickLength(0.05)
    rp_th1 .GetXaxis().SetLabelSize( 0.07)
    rp_th1 .GetYaxis().SetLabelSize( 0.06)

    rp_th1.GetYaxis().SetRangeUser(0.0,3.2)
# #     rp_th1.SetMaximum(2.1)
# #     rp_th1.SetMinimum(0.05)

    lowerPad.cd()
    rp_th1. Draw("")
    line =ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(ROOT.kGreen+4)
    line.Draw()
    r_b .Draw('same P')

    upperPad.cd()
    txt = ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.042)
    txt.DrawLatexNDC(.73, .91, '2025 (13.6 TeV)')   

    c1.Update()
    c1.Modified()

    c1.SaveAs('compare_for_nano/%s.pdf' %(iplot.label))
    c1.SaveAs('compare_for_nano/%s.png' %(iplot.label))
