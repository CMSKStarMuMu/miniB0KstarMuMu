import ROOT
import sys
import math
from DataFormats.FWLite import Events
from collections import OrderedDict
from array import array


b0_mass    = 5.27963
b0_width   = 0.03601 
kaonMass_  = 0.493677
kstMass_   = 0.896

mc_sigma   = 0.0490
mc_mass    = 5.27783 

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

class PlotContainer(object):
    def __init__(self, var, selection, mcselection, xtitle, norm, nbins, xmin, xmax, pltopt, label=None, logx=False, logy=False, fill=False):
        self.var         = var 
        self.selection   = selection 
        self.mcselection = mcselection 
        self.xtitle      = xtitle 
        self.norm        = norm 
        self.nbins       = nbins 
        self.xmin        = xmin 
        self.xmax        = xmax
        self.pltopt      = pltopt
        self.logy        = logy
        self.logx        = logx
        self.fill        = fill
        if label:
            self.label = label
        else:
            self.label = self.var

## data
t1 = ROOT.TChain('ntuple')
t1.Add('../ntuples/out_flat_LMNR_2016H.root')


## mc 
t2 = ROOT.TChain('ntuple')
t2.Add('../ntuples/mc/out_flat_MC_sub1_sub2_sub3_fix_addW.root')


selectionMC   = '((tagB0==1 && (bMass    > {M}-2.5*{S} && bMass    < {M}+2.5*{S}) ) || \
                  (tagB0==0 && (bBarMass > {M}-2.5*{S} && bBarMass < {M}+2.5*{S}) )) && \
                  (mumNTrkLayers >= 6        && mupNTrkLayers >= 6)  && \
                  (mumNPixLayers >= 1        && mupNPixLayers >= 1)  && \
                  (mumdxyVtx < 0.3           && mupdxyVtx < 0.3 )    && \
                  (mumdzVtx < 20             && mupdzVtx  < 20  )    && \
                  (mumHighPurity == 1        && mupHighPurity == 1 ) && \
                  (mumTMOneStationTight == 1 && mupTMOneStationTight == 1 ) && \
                  (kstTrkmHighPurity == 1    && kstTrkpHighPurity == 1) && \
                  (kkMass > 1.035) && \
                  !(kstTrkmGlobalMuon==1 && kstTrkmNTrkLayers > 5 && kstTrkmNPixHits > 0) && \
                  !(kstTrkpGlobalMuon==1 && kstTrkpNTrkLayers > 5 && kstTrkpNPixHits > 0) && \
                    truthMatchMum == 1 && truthMatchMup == 1 && truthMatchTrkm == 1 && truthMatchTrkp == 1'.format(M=mc_mass,S=mc_sigma)

selectionData = '((tagB0==1 && ((bMass    > {M}-7.5*{S} && bMass    < {M}-5*{S}) || (bMass    > {M}+5*{S} && bMass    < {M}+7.5*{S}))) || \
                  (tagB0==0 && ((bBarMass > {M}-7.5*{S} && bBarMass < {M}-5*{S}) || (bBarMass > {M}+5*{S} && bBarMass < {M}+7.5*{S}))) ) && \
                  (mumNTrkLayers >= 6        && mupNTrkLayers >= 6)  && \
                  (mumNPixLayers >= 1        && mupNPixLayers >= 1)  && \
                  (mumdxyVtx < 0.3           && mupdxyVtx < 0.3 )    && \
                  (mumdzVtx < 20             && mupdzVtx  < 20  )    && \
                  (mumHighPurity == 1        && mupHighPurity == 1 ) && \
                  (mumTMOneStationTight == 1 && mupTMOneStationTight == 1 ) && \
                  (kstTrkmHighPurity == 1    && kstTrkpHighPurity == 1) && \
                  (kkMass > 1.035) && \
                  !(kstTrkmGlobalMuon==1 && kstTrkmNTrkLayers > 5 && kstTrkmNPixHits > 0) && \
                  !(kstTrkpGlobalMuon==1 && kstTrkpNTrkLayers > 5 && kstTrkpNPixHits > 0)'.format(M=mc_mass,S=mc_sigma)


thebmass   = 'bMass*tagB0 + bBarMass*(1-tagB0)'
thekstmass = 'kstMass*tagB0 + kstBarMass*(1-tagB0)'


toplot = [
#     PlotContainer(
#         var         = 'abs(cand_refit_mez_2-cand_refit_mez_1)', 
#         selection   = '1',
#         mcselection = 'weight', 
#         xtitle      = '|p_{Z,1}^{miss} - p_{Z,2}^{miss}| [GeV]', 
#         norm        = True, 
#         nbins       = 50, 
#         xmin        = 0, 
#         xmax        = 2500,
#         pltopt      = 'HIST',
#         label       = 'pzdiff'
#     ),

    PlotContainer(var = 'bVtxCL'                         , selection = selectionData, mcselection = selectionMC, xtitle = 'B^{0}* vtx CL'                   , norm = True, nbins = 100, xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstVtxCL'                       , selection = selectionData, mcselection = selectionMC, xtitle = 'K* vtx CL'                       , norm = True, nbins = 100, xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuVtxCL'                      , selection = selectionData, mcselection = selectionMC, xtitle = 'dimuon vtx CL'                   , norm = True, nbins = 100, xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kkMass'                         , selection = selectionData, mcselection = selectionMC, xtitle = 'K^{+}K^{-} mass [GeV]'           , norm = True, nbins = 200, xmin =  0  , xmax =  5    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuMass'                       , selection = selectionData, mcselection = selectionMC, xtitle = 'B0 mass [GeV]'                   , norm = True, nbins =  50, xmin =  4  , xmax =  7    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuCosAlphaBS'                 , selection = selectionData, mcselection = selectionMC, xtitle = 'cos#alpha_{BS}'                  , norm = True, nbins = 100, xmin = -1  , xmax =  1    , fill = True, pltopt = 'HIST', logy = True),
    PlotContainer(var = 'mumuLBS'                        , selection = selectionData, mcselection = selectionMC, xtitle = 'L_{BS}^{J/#psi} [cm]'            , norm = True, nbins = 500, xmin =  0  , xmax =  2    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuLBS/mumuLBSE'               , selection = selectionData, mcselection = selectionMC, xtitle = 'L_{BS}^{J/#psi}/#sigma'          , norm = True, nbins = 500, xmin =  0  , xmax =  100  , fill = True, pltopt = 'HIST', label = 'LSigmaJpsi'),
    PlotContainer(var = 'kstTrkpDCABS/kstTrkpDCABSE'     , selection = selectionData, mcselection = selectionMC, xtitle = 'trk^{+} DCA_{BS} significance'   , norm = True, nbins = 100, xmin =  0  , xmax =  2    , fill = True, pltopt = 'HIST', logy = True, label='trkpDCAsignificance'),
    PlotContainer(var = 'kstTrkmDCABS/kstTrkmDCABSE'     , selection = selectionData, mcselection = selectionMC, xtitle = 'trk^{-} DCA_{BS} significance'   , norm = True, nbins = 100, xmin =  0  , xmax =  2    , fill = True, pltopt = 'HIST', logy = True, label='trkmDCAsignificance'),
    PlotContainer(var = 'kstTrkmCL'                      , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkmCL'                       , norm = True, nbins = 50,  xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpCL'                      , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkpCL'                       , norm = True, nbins = 50,  xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmNormChi2'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkmNormChi2'                 , norm = True, nbins = 50,  xmin =  0  , xmax =  30   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpNormChi2'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkpNormChi2'                 , norm = True, nbins = 50,  xmin =  0  , xmax =  30   , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmFracHits'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkmFracHits'                 , norm = True, nbins = 50,  xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpFracHits'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkpFracHits'                 , norm = True, nbins = 50,  xmin =  0  , xmax =  1    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmNPixHits'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkmNPixHits'                 , norm = True, nbins = 20,  xmin =  0  , xmax = 20    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpNPixHits'                , selection = selectionData, mcselection = selectionMC, xtitle = 'kstTrkpNPixHits'                 , norm = True, nbins = 20,  xmin =  0  , xmax = 20    , fill = True, pltopt = 'HIST'),
         
    PlotContainer(var = 'bPt'                            , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(B^{0}) [GeV]'              , norm = True, nbins = 300, xmin =  0  , xmax = 60    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstPt'                          , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(K*) [GeV]'                 , norm = True, nbins = 300, xmin =  0  , xmax = 30    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuPt'                         , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(#mu^{+}#mu^{-}) [GeV]'     , norm = True, nbins = 300, xmin =  0  , xmax = 60    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumPt'                          , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(#mu^{-}) [GeV]'            , norm = True, nbins = 300, xmin =  0  , xmax = 40    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mupPt'                          , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(#mu^{+}) [GeV]'            , norm = True, nbins = 300, xmin =  0  , xmax = 40    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmPt'                      , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(trk^{-}) [GeV]'            , norm = True, nbins = 300, xmin =  0  , xmax = 30    , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpPt'                      , selection = selectionData, mcselection = selectionMC, xtitle = 'p_{T}(trk^{+}) [GeV]'            , norm = True, nbins = 300, xmin =  0  , xmax = 30    , fill = True, pltopt = 'HIST'),
         
    PlotContainer(var = 'bPhi'                           , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(B^{0})'                     , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstPhi'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(K*)'                        , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuPhi'                        , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(#mu^{+}#mu^{-})'            , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumPhi'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(#mu^{-})'                   , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mupPhi'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(#mu^{+})'                   , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmPhi'                     , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(trk^{-})'                   , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpPhi'                     , selection = selectionData, mcselection = selectionMC, xtitle = '#phi(trk^{+})'                   , norm = True, nbins = 128, xmin = -3.2 , xmax =  3.2  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'bEta'                           , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(B^{0})'                     , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstEta'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(K*)'                        , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumuEta'                        , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(#mu^{+}#mu^{-})'            , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mumEta'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(#mu^{-})'                   , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'mupEta'                         , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(#mu^{+})'                   , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkmEta'                     , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(trk^{-})'                   , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
    PlotContainer(var = 'kstTrkpEta'                     , selection = selectionData, mcselection = selectionMC, xtitle = '#eta(trk^{+})'                   , norm = True, nbins = 200, xmin = -2.5 , xmax =  2.5  , fill = True, pltopt = 'HIST'),
]

c1 = ROOT.TCanvas('c1', 'c1', 700, 700)

for i, iplot in enumerate(toplot):
    
    nbins  = iplot.nbins
    xmin   = iplot.xmin  
    xmax   = iplot.xmax  
    
    sel    = iplot.selection
    mcsel  = iplot.mcselection

    var    = iplot.var
    
    h_data = ROOT.TH1F('h_data_%d'%i, '', nbins, xmin, xmax)
    h_mc   = ROOT.TH1F('h_mc_%d'  %i, '', nbins, xmin, xmax)

    h_data.SetLineColor(ROOT.kBlack)
    h_mc  .SetLineColor(ROOT.kCyan+2)

    xtitle = iplot.xtitle

    h_data.GetXaxis().SetTitle(xtitle)
    h_mc  .GetXaxis().SetTitle(xtitle)

    h_data.GetYaxis().SetTitle('events')
    h_mc  .GetYaxis().SetTitle('events')

    h_data.GetXaxis().SetTitleOffset(1.2)
    h_mc  .GetXaxis().SetTitleOffset(1.2)

    h_data.GetYaxis().SetTitleOffset(1.5)
    h_mc  .GetYaxis().SetTitleOffset(1.5)
    
    t1.Draw('%s >> %s' %(var, h_data.GetName()), '%s'          %(sel), iplot.pltopt         )
    t2.Draw('%s >> %s' %(var, h_mc  .GetName()), '(%s)*weight' %( mcsel), iplot.pltopt + 'SAME')
#     t2.Draw('%s >> %s' %(var, h_mc  .GetName()), '(%s)*(%s)' %(sel, mcsel), iplot.pltopt + 'SAME')

    if iplot.norm:
        h_data.Scale(1./h_data.Integral())
        h_mc  .Scale(1./h_mc  .Integral())

        h_data.GetYaxis().SetTitle('a.u.')
        h_mc  .GetYaxis().SetTitle('a.u.')

    ymax = max(h_data.GetMaximum(), h_mc.GetMaximum())

    h_data.SetMaximum(1.1*ymax)
    h_mc  .SetMaximum(1.1*ymax)
            
    h_data.Draw(iplot.pltopt)
    h_mc  .Draw(iplot.pltopt + 'SAME')

    l1 = ROOT.TLegend(0.56,0.78,0.84,0.88)
    l1.AddEntry(h_data, 'background'  , 'f')
    l1.AddEntry(h_mc  , 'signal'      , 'f')
    l1.SetBorderSize(0)
    l1.Draw('same')

    if iplot.fill:
        h_data.SetFillStyle(0)
        h_mc  .SetFillStyle(3001)
        h_data.SetFillColor(ROOT.kBlack)
        h_mc  .SetFillColor(ROOT.kCyan+2  )
        h_data.SetLineWidth(2)

    ROOT.gPad.SetLogx(iplot.logx)
    ROOT.gPad.SetLogy(iplot.logy)
    
    ROOT.gPad.Update()
    
    c1.SaveAs('data_mc_comparison/%s.pdf' %iplot.label)


