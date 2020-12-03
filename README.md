# miniB0KstarMuMu

**ntuple production**  
```
cmsrel CMSSW_9_4_0_patch1 
cd CMSSW_9_4_0_patch1/src

cmsenv  
git cms-addpkg PhysicsTools/PatAlgos
git remote add packed-cmssw git@github.com:sarafiorendi/cmssw.git
git fetch packed-cmssw
git checkout packed-cmssw/addPATSelector93X PhysicsTools/PatAlgos/plugins/PATObjectSelector.cc
git checkout packed-cmssw/addPATSelector93X PhysicsTools/PatAlgos/plugins/PATObjectSelector.h

git clone git@github.com:CMSKStarMuMu/miniB0KstarMuMu.git 
scram b

cd miniB0KstarMuMu/miniKstarMuMu/test
cmsRun miniKstarMuMu_cfg.py
```

For 10X releases you should instead use
```
git checkout packed-cmssw/addPATSelector1012 PhysicsTools/PatAlgos/plugins/PATObjectSelector.cc
git checkout packed-cmssw/addPATSelector1012 PhysicsTools/PatAlgos/plugins/PATObjectSelector.h
```
**test/flat_ntuples**  
scripts to convert ntuples with array as branches into plain ntuples; multiple candidates per event are kept and converted into multiple entries

**test/flat_ntuples/plot**  
scripts for plotting distributions starting from flat ntuples  

