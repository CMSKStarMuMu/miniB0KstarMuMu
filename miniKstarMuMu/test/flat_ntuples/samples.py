# sample = {}

genericPath = '/gwteras/cms/store/user/fiorendi/p5prime/'

LMNR_sub_A  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018C/190430_081912/'  ## wrong name during crab submission
LMNR_sub_B  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018B/190430_081918/'
LMNR_sub_C  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018A/190430_081926/'  ## wrong name during crab submission
LMNR_sub_D  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018D/190430_081940/'


PSI_sub_A  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018A/190430_090004/'
PSI_sub_B  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018B/190430_090019/'
PSI_sub_C  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018C/190430_090013/'
PSI_sub_D  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018D/190430_085958/'

# 
# samples['LMNR_2018A' ] = [ genericPath + LMNR_sub_A ]
# samples['LMNR_2018B' ] = [ genericPath + LMNR_sub_B ]
# samples['LMNR_2018C' ] = [ genericPath + LMNR_sub_C ]
# samples['LMNR_2018D' ] = [ genericPath + LMNR_sub_D ]
# samples
# samples['PSI_2018A' ] = [ genericPath + PSI_sub_A ]
# samples['PSI_2018B' ] = [ genericPath + PSI_sub_B ]
# samples['PSI_2018C' ] = [ genericPath + PSI_sub_C ]
# samples['PSI_2018D' ] = [ genericPath + PSI_sub_D ]

## ------------------------------
# era_dict['MC_LMNR' ] = [ genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/miniB0KstarMuMu_2018MC_9Mar20_KstWindow0p250_MC_LMNR_Autumn18/200309_130844/' ]

samples = {
    ## Data
    'LMNR_2018A' : {
        'path' : genericPath + LMNR_sub_A,
        'isMC' : False,
    },
    'LMNR_2018B' : {
        'path' : genericPath + LMNR_sub_B,
        'isMC' : False,
    },

    'LMNR_2018C' : {
        'path' : genericPath + LMNR_sub_C,
        'isMC' : False,
    },
    'LMNR_2018D' : {
        'path' : genericPath + LMNR_sub_D,
        'isMC' : False,
    },


    ## Data
    'Charmonium_2018A' : {
        'path' : genericPath + PSI_sub_A,
        'isMC' : False,
    },
    'Charmonium_2018B' : {
        'path' : genericPath + PSI_sub_B,
        'isMC' : False,
    },

    'Charmonium_2018C' : {
        'path' : genericPath + PSI_sub_C,
        'isMC' : False,
    },
    'Charmonium_2018D' : {
        'path' : genericPath + PSI_sub_D,
        'isMC' : False,
    },

    ## MC
    'MC_LMNR_2018' : {
        'path' : genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/miniB0KstarMuMu_2018MC_Apr30_MC_LMNR_Autumn18/190502_085626/',
        'isMC' : True,
    },
   'MC_JPSI_2018' : {
        'path' : genericPath + 'BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/miniB0KstarMuMu_2018MC_Apr30_MC_JPSI_Autumn18/190502_085631/',
        'isMC' : True,
    },

    'MC_LMNR_2017' : {
        'path' : genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/miniB0KstarMuMu_Mar12_MuTk_saveTrigObj_FixIso_MC_LMNR_Fall17MiniAODv2/190416_153857/',
        'isMC' : True,
    },

    'MC_LMNR_2016_p1' : {
        'path' : genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_B0KstarMuMu_part1_asFeb5_LS/190511_090224/',
        'isMC' : True,
    },
    'MC_LMNR_2016_p2' : {
        'path' : genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_B0KstarMuMu_part2_asFeb5_LS/190510_100225/',
        'isMC' : True,
    },
    'MC_LMNR_2016_p3' : {
        'path' : genericPath + 'BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_B0KstarMuMu_part3_asFeb5_LS/190511_092617/',
        'isMC' : True,
    },
#    'MC_PSI' : {
#         'path' : genericPath + 'BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/miniB0KstarMuMu_2018MC_Apr30_MC_JPSI_Autumn18/190502_085631/',
#         'isMC' : True,
#     },

    'MC_LMNR_NOFILTER' : {
        'path' : genericPath + 'BdToKstarMuMu_BFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_GEN_BdMuMuKStar_NoFilter/181204_162315/',
        'isMC' : True,
    },


}