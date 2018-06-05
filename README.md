<<<<<<< HEAD
# miniB0KstarMuMu

**ntuple production**  
```
cmsrel CMSSW_9_3_5 
cd CMSSW_9_3_5_patch1/src

git init
git clone git@github.com:sarafiorendi/miniB0KstarMuMu.git .
scram b

cd miniB0KstarMuMu/miniKstarMuMu/test
cmsRun miniKstarMuMu_cfg.py
```

**test/flat_ntuples**  
scripts to convert ntuples with array as branches into plain ntuples; multiple candidates per event are kept and converted into multiple entries

**test/flat_ntuples/plot**  
scripts for plotting distributions starting from flat ntuples  
=======
# miniB0KstarMuMu
>>>>>>> mib/master
