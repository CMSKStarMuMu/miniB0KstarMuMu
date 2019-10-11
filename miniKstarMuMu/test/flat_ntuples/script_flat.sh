#!/bin/tcsh -f
    source /gwpool/initcms/root-standalone.csh
#             setenv KRB5CCNAME /gwpool/users/fiorendi/krb5cc_30006
#             eosfusebind
    set W_DIR = "/gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src/miniB0KstarMuMu/miniKstarMuMu/test/flat_ntuples"
    cd $W_DIR
    eval `scramv1 runtime -csh`
    python flatNtuplesXX_batch.py -n thenumber -f thefolder dogen -e theera -c thechan
#     python flatNtuplesXX_batch.py -n thenumber -f thefolder -d inputdir