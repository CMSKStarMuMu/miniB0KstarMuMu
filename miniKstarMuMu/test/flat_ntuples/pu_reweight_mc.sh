#!/bin/tcsh

echo '==== running step 1 ===='

setenv yfolder  '/gwteray/users/fiorendi/data2018/p5prime/flat_ntuples/Apr30/'
setenv ifolder  'ntuples/mc/all_gen_psi/'
setenv jfolder  'ntuples/mc/mc_psi/'
setenv genmc    'gen_ntuple_2018MC_PSI_AddVtx.root'
setenv recomc   'reco_ntuple_2018MC_PSI.root'


python addNVtxWeight.py \
../ntuple_data/MyDataPileupHistogram_2018_may6.root  \
$ifolder \
$yfolder$genmc \
-t ntuple


hadd $yfolder$recomc $jfolder/reco*.root

echo 'now adding weights to the reco sample'
python addNVtxWeightReco.py \
$yfolder$genmc \
$yfolder$recomc \
-t ntuple


### obtain pileup file from data doing:
## pileupCalc.py -i total_lumi_processed_apr16.json --inputLumiJSON ../pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 1000 --numPileupBins 1000 MyDataPileupHistogram_2017_apr16.root
#####

# python addNVtxWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2017.root  \
# "/gwteray/users/fiorendi/data2017/p5prime/flat_ntuples/Mar12_hltMatch/ntuple_2017MC_LMNR_newHLTMatch.root" \
# /gwteray/users/fiorendi/data2017/p5prime/flat_ntuples/Mar12_hltMatch/ntuple_2017MC_LMNR_addNVtx.root \
# -t ntuple


# python addNVtxWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2017.root  \
# "/gwteray/users/fiorendi/data2017/p5prime/flat_ntuples/Mar12_hltMatch/ntuple_2017MC_JPSI_newHLTMatch.root" \
# /gwteray/users/fiorendi/data2017/p5prime/flat_ntuples/Mar12_hltMatch/ntuple_2017MC_JPSI_newHLTMatch_AddVtx.root \
# -t ntuple

# echo '==== running step 2 ===='
# python addNVtxWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2016B_2016F.root  \
# /gwteray/users/fiorendi/data2016/p5prime/flat_ntuples/Feb5_skimSoftMu/ntuple_MC_LMNR_PostRefitMomenta_new_addMissingFiles_new_addVtx0.root \
# /gwteray/users/fiorendi/data2016/p5prime/flat_ntuples/Feb5_skimSoftMu/ntuple_MC_LMNR_PostRefitMomenta_new_addMissingFiles_new_addVtx1.root \
# -t ntuple\
# --wName weightBF
# 
# echo '==== running step 3 ===='
# 
# python addNVtxWeight.py \
# ../ntuple_data/MyDataPileupHistogram_2016G_2016H.root  \
# /gwteray/users/fiorendi/data2016/p5prime/flat_ntuples/Feb5_skimSoftMu/ntuple_MC_LMNR_PostRefitMomenta_new_addMissingFiles_new_addVtx1.root \
# /gwteray/users/fiorendi/data2016/p5prime/flat_ntuples/Feb5_skimSoftMu/ntuple_MC_LMNR_PostRefitMomenta_new_addMissingFiles_new_addVtx.root \
# -t ntuple\
# --wName weightGH
# 
