# miniB0KstarMuMu

**ntuple production for Run 3**  
```
cmsrel CMSSW_14_2_2 
cd CMSSW_14_2_2/src

cmsenv  

git clone -b run3_validate_bphnano git@github.com:CMSKStarMuMu/miniB0KstarMuMu.git 
scram b -j

cd miniB0KstarMuMu/miniKstarMuMu/test
cmsRun miniKstarMuMu_cfg.py
```

**test/flat_ntuples**  
scripts to convert ntuples with array as branches into plain ntuples; multiple candidates per event are kept and converted into multiple entries

**test/flat_ntuples/plot**  
scripts for plotting distributions starting from flat ntuples  

