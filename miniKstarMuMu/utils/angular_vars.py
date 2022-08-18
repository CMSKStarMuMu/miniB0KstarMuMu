import ROOT, math
import numpy as np
from ROOT import TLorentzVector
from copy import deepcopy as dc
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi


muonmass_ = 0.1056583745
kaonmass_ = 0.493677
pionmass_ = 0.139570


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
mum_lv          = TLorentzVector()
mup_lv          = TLorentzVector()
tkp_lv          = TLorentzVector()
tkm_lv          = TLorentzVector()
mu_lv_boosted   = TLorentzVector()
kaon_lv_boosted = TLorentzVector()
tk1_lv          = TLorentzVector()
tk2_lv          = TLorentzVector()


def addGenVars(genSignal,
           mumPx,      mumPy,     mumPz,  
           mupPx,      mupPy,     mupPz,  
           kstTrkmPx,  kstTrkmPy, kstTrkmPz,
           kstTrkpPx,  kstTrkpPy, kstTrkpPz,
           kstPx,      kstPy,     kstPz,   kstMass,
           bPx,        bPy,       bPz,     bMass
          ):

  tagB0 = 1 if (genSignal==1 or genSignal==3 or genSignal==5) else 0        
  
  kst_e  = math.sqrt( kstPx**2 + kstPy**2 + kstPz**2 + kstMass**2 )
  kst_lv.SetPxPyPzE(kstPx, kstPy, kstPz, kst_e)

  b0_e  = math.sqrt( bPx**2 + bPy**2 + bPz**2 + bMass**2 )
  b0_lv.SetPxPyPzE(bPx, bPy, bPz, b0_e)

  mum_e  = math.sqrt( mumPx**2 + mumPy**2 + mumPz**2 + muonmass_**2 )
  mum_lv.SetPxPyPzE(mumPx, mumPy, mumPz, mum_e)

  mup_e  = math.sqrt( mupPx**2 + mupPy**2 + mupPz**2 + muonmass_**2 )
  mup_lv.SetPxPyPzE(mupPx, mupPy, mupPz, mup_e)

  mumu_lv = mup_lv + mum_lv

  tkp_e  = math.sqrt( kstTrkpPx**2 + kstTrkpPy**2 + kstTrkpPz**2 + tagB0*kaonmass_**2 + (1-tagB0)*pionmass_**2 )
  tkp_lv.SetPxPyPzE(kstTrkpPx, kstTrkpPy, kstTrkpPz, tkp_e)

  tkm_e  = math.sqrt( kstTrkmPx**2 + kstTrkmPy**2 + kstTrkmPz**2 + tagB0*pionmass_**2 + (1-tagB0)*kaonmass_**2 )
  tkm_lv.SetPxPyPzE(kstTrkmPx, kstTrkmPy, kstTrkmPz, tkm_e)

  cos_theta_k     = -99.
  cos_theta_l     = -99.
  phiKstMuMuPlane = -99.
  mu_lv_boosted   = TLorentzVector()
  kaon_lv_boosted = TLorentzVector()

  if tagB0:
    mu_lv_boosted   = dc(mup_lv)
    kaon_lv_boosted = dc(tkp_lv)
  
  else:
    mu_lv_boosted   = dc(mum_lv)
    kaon_lv_boosted = dc(tkm_lv)
    
  ## calculate cos theta_l
  boostMuMu = mumu_lv.BoostVector()
  mu_lv_boosted.Boost(-boostMuMu)
  b0_lv_boosted = dc(b0_lv)
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
  b0_lv_boosted = dc(b0_lv)
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
  mum_lv_boosted = dc(mum_lv)
  mum_lv_boosted.Boost(-boostB0)

  mup_lv_boosted = dc(mup_lv)
  mup_lv_boosted.Boost(-boostB0)
  
  tkp_lv_boosted = dc(tkp_lv)
  tkp_lv_boosted.Boost(-boostB0)

  tkm_lv_boosted = dc(tkm_lv)
  tkm_lv_boosted.Boost(-boostB0)
  
  MuMuPlane = mup_lv_boosted.Vect().Cross(mum_lv_boosted.Vect())
  KstPlane  = tkp_lv_boosted.Vect().Cross(tkm_lv_boosted.Vect())

  kst_lv_boosted = dc(kst_lv)
  kst_lv_boosted.Boost(-boostB0)
  phiKstMuMuPlane = -99

  if tagB0==1:
    if MuMuPlane.Cross(KstPlane).Dot(kst_lv_boosted.Vect()) > 0:
      phiKstMuMuPlane = MuMuPlane.Angle(KstPlane)
    else:
      phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane)
  elif tagB0==0:
    if MuMuPlane.Cross(KstPlane).Dot(kst_lv_boosted.Vect()) > 0:
      phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane)
    else:
      phiKstMuMuPlane = MuMuPlane.Angle(KstPlane)

  return cos_theta_l, cos_theta_k, phiKstMuMuPlane#, mumu_lv.M()
  
  
  
  
def addVars(tagB0,
           mumuPx,     mumuPy,    mumuPz,  mumuMass,
           mumPx,      mumPy,     mumPz,  
           mupPx,      mupPy,     mupPz,  
           kstTrkmPx,  kstTrkmPy, kstTrkmPz,
           kstTrkpPx,  kstTrkpPy, kstTrkpPz,
           kstPx,      kstPy,     kstPz,   kstMass, kstBarMass,
           bPx,        bPy,       bPz,     bMass,   bBarMass
          ):
          
          
  if mumuPx == -99:
      return -99, -99, -99        
                  
  mumu_e  = math.sqrt( mumuPx**2 + mumuPy**2 + mumuPz**2 + mumuMass**2 )
  mumu_lv.SetPxPyPzE(mumuPx, mumuPy, mumuPz, mumu_e)

  kst_e  = math.sqrt( kstPx**2 + kstPy**2 + kstPz**2 + tagB0*kstMass**2 + (1-tagB0)*kstBarMass**2 )
  kst_lv.SetPxPyPzE(kstPx, kstPy, kstPz, kst_e)

  b0_e  = math.sqrt( bPx**2 + bPy**2 + bPz**2 + tagB0*bMass**2 + (1-tagB0)*bBarMass**2 )
  b0_lv.SetPxPyPzE(bPx, bPy, bPz, b0_e)

  mum_e  = math.sqrt( mumPx**2 + mumPy**2 + mumPz**2 + muonmass_**2 )
  mum_lv.SetPxPyPzE(mumPx, mumPy, mumPz, mum_e)

  mup_e  = math.sqrt( mupPx**2 + mupPy**2 + mupPz**2 + muonmass_**2 )
  mup_lv.SetPxPyPzE(mupPx, mupPy, mupPz, mup_e)

  tkp_e  = math.sqrt( kstTrkpPx**2 + kstTrkpPy**2 + kstTrkpPz**2 + tagB0*kaonmass_**2 + (1-tagB0)*pionmass_**2 )
  tkp_lv.SetPxPyPzE(kstTrkpPx, kstTrkpPy, kstTrkpPz, tkp_e)

  tkm_e  = math.sqrt( kstTrkmPx**2 + kstTrkmPy**2 + kstTrkmPz**2 + tagB0*pionmass_**2 + (1-tagB0)*kaonmass_**2 )
  tkm_lv.SetPxPyPzE(kstTrkmPx, kstTrkmPy, kstTrkmPz, tkm_e)
  
  cos_theta_k     = -99.
  cos_theta_l     = -99.
  phiKstMuMuPlane = -99.

  if tagB0:
    mu_lv_boosted   = dc(mup_lv)
    kaon_lv_boosted = dc(tkp_lv)
  
  else:
    mu_lv_boosted = dc(mum_lv)
    kaon_lv_boosted = dc(tkm_lv)
    
  ## calculate cos theta_l
  boostMuMu = mumu_lv.BoostVector()
  mu_lv_boosted.Boost(-boostMuMu)
  b0_lv_boosted = dc(b0_lv)
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
  b0_lv_boosted = dc(b0_lv)
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
  mum_lv_boosted = dc(mum_lv)
  mum_lv_boosted.Boost(-boostB0)

  mup_lv_boosted = dc(mup_lv)
  mup_lv_boosted.Boost(-boostB0)
  
  tkp_lv_boosted = dc(tkp_lv)
  tkp_lv_boosted.Boost(-boostB0)

  tkm_lv_boosted = dc(tkm_lv)
  tkm_lv_boosted.Boost(-boostB0)
  
  MuMuPlane = mup_lv_boosted.Vect().Cross(mum_lv_boosted.Vect())
  KstPlane  = tkp_lv_boosted.Vect().Cross(tkm_lv_boosted.Vect())

  kst_lv_boosted = dc(kst_lv)
  kst_lv_boosted.Boost(-boostB0)
  phiKstMuMuPlane = -99
  if tagB0==1:
    if MuMuPlane.Cross(KstPlane).Dot(kst_lv_boosted.Vect()) > 0:
      phiKstMuMuPlane = MuMuPlane.Angle(KstPlane)
    else:
      phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane)
  elif tagB0==0:
    if MuMuPlane.Cross(KstPlane).Dot(kst_lv_boosted.Vect()) > 0:
      phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane)
    else:
      phiKstMuMuPlane = MuMuPlane.Angle(KstPlane)

  return cos_theta_l, cos_theta_k, phiKstMuMuPlane




def addMMKVars(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkpPt,  tkpEta,  tkpPhi,
            tkmPt,  tkmEta,  tkmPhi
          ):
          
    if mumPt == -99:
        return -99, -99, -99     
    
    mmk1 = -99.
    mmk2 = -99.
        
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
            
    if tkpPt >= tkmPt:                
      tk1_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
      tk2_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
    else:                
      tk1_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
      tk2_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
  
    mmk1 = (mum_lv+mup_lv+tk1_lv).M()
    mmk2 = (mum_lv+mup_lv+tk2_lv).M()
    
    return mmk1, mmk2


def addDR(
           mumEta,      mumPhi,     kstTrkmEta,  kstTrkmPhi
          ):
  return deltaR(mumEta, mumPhi, kstTrkmEta, kstTrkmPhi)        
