#!/bin/tcsh

##################
## Updated usage

## - first subit gen-only ntuples to prepare the trueNI histogram
##   this can be done using the flatNtuplesMCWeight_batch.py script and the usual submission
## - then run this script to prepare the weight histogram
## - finally run your MC ntuples as usual using flatNtuplesMC_batch.py
##################

echo '==== running step 1 ===='

source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

setenv ifolder    'ntuples/mc/calcWeight_2017_bs/'
setenv year        2017
setenv channel     'BS'


######### for 2016 #######################
# python computePileupWeight.py \
# "../ntuple_data/MyDataPileupHistogram_2016B_2016H.root ../ntuple_data/MyDataPileupHistogram_2016B_2016F.root ../ntuple_data/MyDataPileupHistogram_2016G_2016H.root"  \
# "$ifolder" \
# $year \
# -t ntuple \
# -c $channel

## 2017 
python computePileupWeight.py \
../ntuple_data/MyDataPileupHistogram_2017_apr16_800bins.root  \
"$ifolder" \
$year \
-t ntuple \
--histogramRange "800,0,800" \
-c $channel 

##  2018
# python computePileupWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2018_may6.root  \
# "$ifolder" \
# $year \
# -t ntuple \
# -c $channel 


## old version
# python addNVtxWeight2016.py \
# "../ntuple_data/MyDataPileupHistogram_2016B_2016H.root ../ntuple_data/MyDataPileupHistogram_2016B_2016F.root ../ntuple_data/MyDataPileupHistogram_2016G_2016H.root"  \
# "$ifolder" \
# $yfolder$genmc \
# -t ntuple
