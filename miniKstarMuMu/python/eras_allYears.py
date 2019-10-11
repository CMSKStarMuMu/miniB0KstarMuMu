eras      = {}
lumi_eras = {}  ## fb-1
run_eras  = {}
 

eras['2016B'] = [ 272007, 275376 ]
eras['2016C'] = [ 275657, 276283 ]
eras['2016D'] = [ 276315, 276811 ]
eras['2016E'] = [ 276831, 277420 ]
eras['2016F'] = [ 277772, 278808 ]
eras['2016G'] = [ 278820, 280385 ]
eras['2016H'] = [ 280919, 284044 ]

eras['APVL1'] = [ 272007, 277990 ]  
eras['APV'  ] = [ 277991, 278801 ]
eras['GOOD' ] = [ 278802, 284044 ]


lumi_eras['2016B']  =  5.805
lumi_eras['2016C']  =  2.615
lumi_eras['2016D']  =  4.281
lumi_eras['2016E']  =  4.035
lumi_eras['2016F']  =  3.119
lumi_eras['2016G']  =  7.658
lumi_eras['2016H']  =  8.773

lumi_eras['APVL1']  = 16.741
lumi_eras['APV'  ]  =  2.7
lumi_eras['GOOD' ]  = 16.831

lumi_eras['2016']  =  lumi_eras['2016B'] + \
                      lumi_eras['2016C'] + \
                      lumi_eras['2016D'] + \
                      lumi_eras['2016E'] + \
                      lumi_eras['2016F'] + \
                      lumi_eras['2016G'] + \
                      lumi_eras['2016H'] 

run_eras[272007,275376] = '2016B' 
run_eras[275657,276283] = '2016C' 
run_eras[276315,276811] = '2016D' 
run_eras[276831,277420] = '2016E' 
run_eras[277772,278808] = '2016F' 
run_eras[278820,280385] = '2016G' 
run_eras[280919,284044] = '2016H' 

run_eras[272007,277990] = 'APVL1' 
run_eras[277991,278801] = 'APV' 
run_eras[278802,284044] = 'GOOD' 


###################################################
## 2017
###################################################
eras['2017B'] = [ 297046, 299329 ]
eras['2017C'] = [ 299368, 302029 ]
eras['2017D'] = [ 302030, 303434 ]
eras['2017E'] = [ 303824, 304797 ]
eras['2017F'] = [ 305040, 306462 ]

eras['NODCDC'] = [ 299368, 304504 ]  ## starting from 2017C to exclude trigger & commissioning problems
eras['DCDC'  ] = [ 304505, 306462 ]

eras['2017D_normalFill'] = [ 302030, 302317 ]  
eras['2017D_8b4e'      ] = [ 302318, 303434 ]  

lumi_eras['2017B']  =   4.942
lumi_eras['2017C']  =   9.840
lumi_eras['2017D']  =   4.280
lumi_eras['2017E']  =   9.413
lumi_eras['2017F']  =  13.602

lumi_eras['NODCDC'] =  20.035# was  24.977 including 2017B
lumi_eras['DCDC'  ] =  17.099

lumi_eras['2017']  =  lumi_eras['2017B'] + \
                      lumi_eras['2017C'] + \
                      lumi_eras['2017D'] + \
                      lumi_eras['2017E'] + \
                      lumi_eras['2017F'] 


run_eras[297046,299329] = '2017B' 
run_eras[299368,302029] = '2017C' 
run_eras[302030,303434] = '2017D' 
run_eras[303824,304797] = '2017E' 
run_eras[305040,306462] = '2017F' 

run_eras[299368,304504] = 'NODCDC' 
run_eras[304505,306462] = 'DCDC' 

# eras['NODCDC'] = [ 297046, 304504 ]
# eras['DCDC'  ] = [ 304505, 306462 ]



###################################################
## 2018
###################################################

eras['2018A'] = [ 315252, 316995 ]
eras['2018B'] = [ 317080, 319310 ]
eras['2018C'] = [ 319337, 320065 ]
eras['2018D'] = [ 320673, 325175 ]

lumi_eras['2018A']  =   14.741
lumi_eras['2018B']  =   7.149
lumi_eras['2018C']  =   6.899
lumi_eras['2018D']  =   32.347

lumi_eras['2018']  =  lumi_eras['2018A'] + \
                      lumi_eras['2018B'] + \
                      lumi_eras['2018C'] + \
                      lumi_eras['2018D'] 


run_eras[315252,316995] = '2018A' 
run_eras[317080,319310] = '2018B' 
run_eras[319337,320065] = '2018C' 
run_eras[320673,325175] = '2018D' 




###################################################
## MC
###################################################
lumi_mc = {}  ## fb-1

lumi_mc['LMNR2016']  =   6216
lumi_mc['JPSI2016']  =     45.6
lumi_mc['PSI2016' ]  =     52.8 ### FIXME

lumi_mc['LMNR2017']  =   8957
lumi_mc['JPSI2017']  =     51.1
lumi_mc['PSI2017' ]  =     93.5

lumi_mc['LMNR2018']  =   8256
lumi_mc['JPSI2018']  =     41.5
lumi_mc['PSI2018' ]  =     75.6

















