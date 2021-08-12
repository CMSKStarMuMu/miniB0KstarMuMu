# Flat ntuples
Set of scripts to skim and flatten the ntuples created by the miniKstarMuMu package.

## How to run
The location of the input ntuples is collected in 
[samples.py](https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/test/flat_ntuples/samples.py).

The two scripts performing the skimming/flattening are 
[flatNtuples_batch.py](https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/test/flat_ntuples/flatNtuples_batch.py) and [flatNtuplesMC_batch.py](https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/test/flat_ntuples/flatNtuplesMC_batch.py), 
and they can be submitted to the queue system (using condor) via the 
[submit_condor.py](https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/test/flat_ntuples/submit_condor.py) script.

```
python submit_condor.py --help
```
will print the input parameters needed by the script to run.
You will find the flat ntuples in the folder ntuple/outdir, where outdir is one of the input parameters to the submit_condor.py script.

## A note on the PU reweight for the MC

In order to add to the MC ntuples (GEN and RECO) the weight needed to reproduce the pileup distribution in data, the script [pu_reweight_mc.sh](https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/test/flat_ntuples/pu_reweight_mc.sh) is available.
It requires that you have first prepared gen-only ntuples to fill the trueNI histogram.
These ntuples are those created by the flatNtuplesMCWeight\_batch.py script, which can be submitted to condor as explained in the previous paragraph. 
Simply run it with 
```
./pu_reweight_mc.sh
```
after having changed the input/output folders appropriately.

