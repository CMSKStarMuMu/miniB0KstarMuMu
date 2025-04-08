import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.MessageLogger = cms.Service("MessageLogger",
# 
#      destinations =  cms.untracked.vstring("long_job_report", 
#      			                            ),
#      long_job_report = cms.untracked.PSet(                       
#              threshold = cms.untracked.string( 'INFO' )          
#      )
# )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#
# It might be better to set the GlobalTag explicitly in productions, otherwise there will be unpredictable effects due to the 
# automatic selection of the GlobalTag. Since this script is just for testing, we can leave it as is. 
#
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data')
#process.GlobalTag.globaltag = cms.string('123X_dataRun3_HLT_v14') ## for data

#process.GlobalTag.globaltag = cms.string('102X_dataRun2_v14') ## for data

process.source = cms.Source("PoolSource",
#BPark
    fileNames = cms.untracked.vstring(
#   'file://19fa6c32-3757-4312-9bbc-35c6a0b7c6bd.root',
#   'file://13d368c5-3e9d-43ce-93fb-62ec35ff0b53.root',
#   'file://010f5c73-45d9-44e9-a1dc-f77f10e09410.root',

#

     'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass1/MINIAOD/10Dec2022-v3/2540000/19fa6c32-3757-4312-9bbc-35c6a0b7c6bd.root',
     'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass1/MINIAOD/10Dec2022-v3/2540000/13d368c5-3e9d-43ce-93fb-62ec35ff0b53.root',
     'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass1/MINIAOD/10Dec2022-v3/2540000/010f5c73-45d9-44e9-a1dc-f77f10e09410.root',
    ),
#     lumisToProcess = cms.untracked.VLuminosityBlockRange('275282:93-275282:94'),
#     eventsToProcess = cms.untracked.VEventRange(
#                                                 '276315:34493463',
#                                                 '276361:515692537',
#                                                 '276776:2211826216',
#     ),
)

taskB0 = cms.Task()

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
taskB0.add(process.unpackedTracksAndVertices)

# process.selectedPFCandidatesHP = cms.EDFilter("PATPackedCandidateSelector",
#     src = cms.InputTag("packedPFCandidates"),
#     cut = cms.string("pt() > 0.8 && abs(eta()) < 2.5 && trackHighPurity() > 0")
# )
# taskB0.add(process.selectedPFCandidatesHP)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')


process.B0KstMuMu = cms.EDAnalyzer("miniKstarMuMu",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons    = cms.InputTag("slimmedMuons"),
    tracks    = cms.InputTag("unpackedTracksAndVertices"),#"selectedPFCandidatesHP"),
#     tracks   = cms.InputTag("selectedPFCandidatesHP"),
    beamSpot = cms.InputTag("offlineBeamSpot"),

    bits          = cms.InputTag("TriggerResults","","HLT"),
    prescales     = cms.InputTag("patTrigger"),
    objects       = cms.InputTag("slimmedPatTrigger"),
    stageL1Trigger = cms.uint32(1),
#    TriggerNames  = cms.vstring("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v", 
#                                "HLT_DoubleMu4_JpsiTrk_Displaced_v",
#                                "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"
    TriggerNames  = cms.vstring("HLT_DoubleMu4_LowMass_Displaced_v", 
                                "HLT_DoubleMu4_Jpsi_Displaced_v",
                                "HLT_DoubleMu4_MuMuTrk_Displaced_v"
                    ),
    L1Names  = cms.vstring("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 
			  "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
			  "L1_DoubleMu4_SQ_OS_dR_Max1p2", 
			  "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
		   ),
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles"),

    PuInfoTag	     = cms.InputTag("slimmedAddPileupInfo"),
    l1results	     = cms.InputTag("gtStage2Digis"),

    ## HLT selections
    MuMuVtxCL	     = cms.untracked.double(0.1 ),  # mu-mu Vtx CL [0.1]
    MuMuLsBS	     = cms.untracked.double(3	),    # mu-mu L/sigma w/respect to BS [3.0]
    DCAMuMu	     = cms.untracked.double(0.5 ),    # mu-mu DCA w/respect to each other [0.5 cm]
    DCAMuBS	     = cms.untracked.double(2.0 ),    # mu DCA w/respect to BS [2.0 cm]
    cosAlphaMuMuBS   = cms.untracked.double(0.9 ),    # mu-mu cos(alpha) w/respect to BS [0.9]
    MinMupT	     = cms.untracked.double(4.0 ),    # mu min pT [4.0 GeV/c]
    MuEta	     = cms.untracked.double(2.4 ),    # mu max eta [2.4]
    MuMupT	     = cms.untracked.double(6.9 ),    # mu-mu min pT [6.9 GeV/c]
    MinMuMuMass      = cms.untracked.double(1.0 ),    # mu-mu min inv. mass [1.0 GeV/c2]
    MaxMuMuMass      = cms.untracked.double(4.8 ),    # mu-mu max inv. mass [4.8 GeV/c2]
    ## Cand pre-selections
     MinB0Mass	     = cms.untracked.double(4.5 ),    # B0 mass lower limit [4.5 GeV/c2]
     MaxB0Mass	     = cms.untracked.double(6.5 ),    # B0 mass upper limit [6.5 GeV/c2]
    B0VtxCL	     = cms.untracked.double(0.01),  # B0 Vtx CL [0.01]
    KstMass	     = cms.untracked.double(3.0 ),    # K*0 (OR K*0bar) mass window sigma [3.0]
    HadDCASBS	     = cms.untracked.double(0.8 ),    # hadron DCA/sigma w/respect to BS [0.8] (also in HLT, now is 2)
    HadpT	     = cms.untracked.double( .8 ),    # hadron min pT [0.8 GeV/c] (also in HLT)
    MaxB0RoughMass   = cms.untracked.double(20. ),    # B0 mass upper limit  before performing the fit#     electrons = cms.InputTag("slimmedElectrons"),

    printMsg	     = cms.untracked.bool(False)
)

process.TFileService = cms.Service('TFileService', 
    fileName = cms.string(
	'B0ToKstMuMu_miniaod.root'
    ), 
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.unpackedTracksAndVertices * process.B0KstMuMu)
# process.p = cms.Path(process.selectedPFCandidatesHP * process.B0KstMuMu)
# process.Tracer = cms.Service("Tracer")

'''
[fiorendi@hercules test]$ dasgoclient -query="parent file=/store/mc/RunIIFall17MiniAODv2/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_N1_94X_mc2017_realistic_v14-v1/20000/08D2B3E9-0043-E811-B981-C4346BC78D10.root"
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/00000/008D78A0-8800-E811-A370-008CFAC91CAC.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/00000/90016980-7C00-E811-A888-008CFAC93EB4.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/00000/9619DB05-5800-E811-94D2-008CFAE45134.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/10000/1C1CD1E3-35F3-E711-9721-E0071B7A8590.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/10000/8ACCBFDE-0FF3-E711-9482-5065F3812271.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/10000/A46C4FE1-0EF3-E711-B34C-24BE05C6B701.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/60000/309C4C23-63F4-E711-A8B2-44A84225CABC.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/60000/9E068997-98F4-E711-B072-44A84225CDA4.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/70001/206A7A4A-65FD-E711-8EA2-FA163E0448E5.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/70001/30E32997-53FD-E711-93F4-FA163EC5F5DC.root
/store/mc/RunIIFall17DRPremix/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/N1_94X_mc2017_realistic_v10-v1/70001/6679548B-54FD-E711-9459-FA163EA32224.root
'''
