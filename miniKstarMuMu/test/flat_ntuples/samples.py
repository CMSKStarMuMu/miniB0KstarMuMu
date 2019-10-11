era_dict = {}

genericPath = '/gwteras/cms/store/user/fiorendi/p5prime/'

LMNR_sub_A  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018C/190430_081912/'  ## wrong name during crab submission
LMNR_sub_B  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018B/190430_081918/'
LMNR_sub_C  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018A/190430_081926/'  ## wrong name during crab submission
LMNR_sub_D  = 'DoubleMuonLowMass/miniB0KstarMuMu_Apr30_LMNR_17Sep2018_2018D/190430_081940/'


PSI_sub_A  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018A/190430_090004/'
PSI_sub_B  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018B/190430_090019/'
PSI_sub_C  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018C/190430_090013/'
PSI_sub_D  = 'Charmonium/miniB0KstarMuMu_Apr30_Charmonium_17Sep2018_2018D/190430_085958/'



era_dict['LMNR_2018A' ] = [ genericPath + LMNR_sub_A ]
era_dict['LMNR_2018B' ] = [ genericPath + LMNR_sub_B ]
era_dict['LMNR_2018C' ] = [ genericPath + LMNR_sub_C ]
era_dict['LMNR_2018D' ] = [ genericPath + LMNR_sub_D ]

era_dict['PSI_2018A' ] = [ genericPath + PSI_sub_A ]
era_dict['PSI_2018B' ] = [ genericPath + PSI_sub_B ]
era_dict['PSI_2018C' ] = [ genericPath + PSI_sub_C ]
era_dict['PSI_2018D' ] = [ genericPath + PSI_sub_D ]

