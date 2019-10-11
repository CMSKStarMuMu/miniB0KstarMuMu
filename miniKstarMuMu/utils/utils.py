muonMass  = 0.10565837
b0_mass   = 5.27963
b0_width  = 0.03601 
kaonMass  = 0.493677

kstMass_  = 0.896
kstSigma_ = 0.05

import ROOT
import math
import numpy as np
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

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


def computeEta(px,py,pz):
    P = math.sqrt(px*px + py*py + pz*pz)
    if ((P - pz) == 0):  
        return -99
    return 0.5*math.log((P + pz) / (P - pz))
    

def computePhi(px,py):
    if px == 0:  
        return -99
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


def findFiringL1(eventL1, eventPrescales, theL1):

    l1list = list(eventL1)
    try:
        ind = l1list.index(theL1)
        pre = eventPrescales[ind]
    except:
        ind = -1
        pre = -1
    return pre


def findTriggerMatching(triplets, 
                        rawmumEta,     rawmumPhi,      rawmumPt,
                        rawmupEta,     rawmupPhi,      rawmupPt,
                        rawkstTrkmEta, rawkstTrkmPhi,  rawkstTrkmPt,
                        rawkstTrkpEta, rawkstTrkpPhi,  rawkstTrkpPt
                       ):

    offline_online_dr = []
    charge_match = 0
    for itr in triplets:
            
        dr_mum  = deltaR(rawmumEta,     rawmumPhi,     itr[0].eta, itr[0].phi)
        dr_mup  = deltaR(rawmupEta,     rawmupPhi,     itr[1].eta, itr[1].phi)
        dr_trkm = deltaR(rawkstTrkmEta, rawkstTrkmPhi, itr[2].eta, itr[2].phi) 
        dr_trkp = deltaR(rawkstTrkpEta, rawkstTrkpPhi, itr[2].eta, itr[2].phi) 
        match_tk_charge = -1      if dr_trkm < dr_trkp else 1
        dr_trk          = dr_trkm if dr_trkm < dr_trkp else dr_trkp

        offline_online_dr.append([dr_mum, dr_mup, dr_trk*match_tk_charge])
        
    ## choose online combination closer in dR, if each dR < 0.1
    best_match    = 0.3
    triplet_index = -1
    for im in range(len(offline_online_dr)):
        if any(abs(t) < 0.1 for t in offline_online_dr[im]):
            sum_dr =  sum_dr_(offline_online_dr[im])
            if sum_dr < best_match:
                best_match   = sum_dr
                charge_match = np.sign(offline_online_dr[im][2])
                triplet_index = im
        
    if triplet_index < 0:
        return charge_match

    return charge_match

sum_dr_ = lambda y : np.sum(np.abs(y))



def muonIso(muiso, num=20):
    '''calculate # iso tracks below threshold'''

    counts = np.zeros(num)   
    if len(muiso) > 0:
        for i in range(num):
            counts[i] = np.size(np.where(np.array(muiso) < 0.1*(i+1)))

    toreturn = ROOT.std.vector('double')() 
    for i in counts:  toreturn.push_back(i)
    return toreturn



