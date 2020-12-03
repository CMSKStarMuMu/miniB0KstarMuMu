#!/bin/tcsh

echo '==== running step 1 ===='

source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

setenv genfolder  '/gwteray/users/fiorendi/p5prime/data2016/flat_ntuples/Oct2020/'
setenv yfolder    '/gwteray/users/fiorendi/p5prime/data2016/flat_ntuples/Oct2020/'
setenv ifolder    'ntuples/mc/mc2016_p1_lumi/'
setenv iifolder   'ntuples/mc/mc2016_p2_lumi/'
setenv iiifolder  'ntuples/mc/mc2016_p3_lumi/'
setenv jfolder    'ntuples/mc/newphi_MC2016_p1/'
setenv genmc      '2016GEN_MC_LMNR.root'
setenv recomc     'reco_ntuple_2016MC_newphi.root'



python addNVtxWeight.py \
../ntuple_data/MyDataPileupHistogram_2018_may6.root  \
$ifolder \
$yfolder$genmc \
-t ntuple

# python addNVtxWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2017_apr16_800bins.root  \
# $ifolder \
# $yfolder$genmc \
# -t ntuple

echo "hadd $yfolder$recomc $jfolder/reco*.root"
hadd $yfolder$recomc $jfolder/reco*.root

echo 'now adding weights to the reco sample'
python addNVtxWeightReco.py \
$genfolder$genmc \
$yfolder$recomc \
-t ntuple

######### for 2016 #######################
# python addNVtxWeight2016.py \
# "../ntuple_data/MyDataPileupHistogram_2016B_2016H.root ../ntuple_data/MyDataPileupHistogram_2016B_2016F.root ../ntuple_data/MyDataPileupHistogram_2016G_2016H.root"  \
# "$ifolder $iifolder $iiifolder" \
# $yfolder$genmc \
# -t ntuple

# echo "hadd $yfolder$recomc $ifolder/reco*.root $iifolder/reco*.root $iiifolder/reco*.root"
# hadd $yfolder$recomc $ifolder/reco*.root $iifolder/reco*.root $iiifolder/reco*.root

# echo 'now adding weights to the reco sample'
# python addNVtxWeightReco.py \
# $genfolder$genmc \
# $yfolder$recomc \
# -t ntuple


### obtain pileup file from data doing:
## pileupCalc.py -i total_lumi_processed_apr16.json --inputLumiJSON ../pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 1000 --numPileupBins 1000 MyDataPileupHistogram_2017_apr16.root
#####

