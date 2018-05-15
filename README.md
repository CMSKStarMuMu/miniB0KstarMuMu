# miniB0KstarMuMu

here using 940patch1, correct release for PhaseII MC is 9_3_5

```
cmsrel CMSSW_9_4_0_patch1 
cd CMSSW_9_4_0_patch1/src

git init
git clone git@github.com:sarafiorendi/miniB0KstarMuMu.git .
scram b

cd miniB0KstarMuMu/miniKstarMuMu/test
cmsRun miniKstarMuMu_cfg.py
