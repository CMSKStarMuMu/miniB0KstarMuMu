# ! /bin/env python

import os
import subprocess
import datetime
from argparse import ArgumentParser
import pdb
import math
from samples import samples
from glob import glob
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument("analyzer", help = "which analyser to run", default = 'BJpsiK_ee_mc' )
parser.add_argument("samples", help = "samples", nargs = '+', choices = samples.keys(), default = 'BJpsiK_ee_mc_2019Oct25' )
parser.add_argument("-n"  , "--njobs"  , dest = "njobs"  , type = int, help = "tot number of input files to be read. All = -1" , default = -1                            )
parser.add_argument("-d"  , "--outdir" , dest = "outdir" ,  help = "output dir"                                     , default = "ntuples" )
parser.add_argument("-a"  , "--addtag" , dest = "addtag" ,  help = "add tag to output dir"                          , default = "ntuples" )
parser.add_argument("-t"  , "--test"   , dest = "test"   ,  help = "do not submit to queue"                        , default = False, action='store_true')
parser.add_argument("--print"          , dest = "printN" ,  help = "print infos"                                   , default = False, action='store_true'      )
parser.add_argument("-S"  , "--start"  , dest = "start"  ,  help = "choose starting file"                           , default =  0                            )
parser.add_argument("-c"  , "--chan"   , dest = "channel",  help = "LMNR, Psi"                                  , default = "LMNR"                 )
parser.add_argument("-g"  , "--dogen"  , dest = "dogen"  ,  help = "produce gen ntuples "                       , default=False, action='store_true'      )
parser.add_argument("-f"  , "--flavour",                    help = "job flavour (sets running time) https://indico.cern.ch/event/731021/contributions/3013463/attachments/1656036/2651022/18-05-24_HTCondor_CMG.pdf", default =  'microcentury', choices = ['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow', 'testmatch', 'nextweek'])
args = parser.parse_args()


script_loc = os.path.realpath(args.analyzer)
channel = args.channel
gen_flag = '--dogen' if (args.dogen) else ''
isMC    = samples[args.samples[0]]['isMC']
mc_flag = '--mc' if isMC==False else ''

base_out = '%s/%s/%s' %(args.outdir, isMC*'mc' + (1-isMC)*'', args.addtag)
# for sample_name in args.samples:
#     sample = samples[sample_name]
    ## local out folder for logs
#     base_out = f'{args.outdir}/{channel}/{sample_name}_{args.addtag}'
os.makedirs('%s/scripts'%base_out)
os.makedirs('%s/outCondor'%base_out)
os.system('cp FindValueFromVectorOfBool.h {base_out}'.format(base_out=base_out))

##  output folder for root files
full_eos_out = '{eos_out_folder}/{base_out}/'.format(eos_out_folder = os.getcwd(), base_out = base_out)
    

## calculate n jobs to be submitted 
import fnmatch
njobs = int(args.njobs   )
if args.njobs == -1:
#     set_trace()
    flist = glob(samples[args.samples[0]]['path']+'/000*/*.root')
    njobs = len(flist)    
#     else:
#         njobs = len(fnmatch.filter(os.listdir(samples['MC_'+args.channel]['path']), '*.root'))    
    
## find missing files from crab
job_n = []
for i,f in enumerate(flist):  job_n.append(int(f.split('_')[-1].split('.')[-2]))    
job_n.sort()    

missing =  sorted(set(range(job_n[0], job_n[-1])) - set(job_n)) 
print 'missing files from crab:', missing
njobs = njobs + len(missing) + 1 ## adding 1 since condor start from 0 and crab_0.root does not exists
print 'n jobs to be submitted: ', njobs 

bname = os.path.realpath('%s/scripts/script_flat.sh'%base_out)


## is nstart!=0 -> instead of process ID maybe use X+processID?
    ## write script_flat.sh script
with open(bname, 'w') as batch:
    batch.write('''#!/bin/tcsh
setenv CMSSWDIR /gwpool/users/fiorendi/p5prime/miniAOD/CMSSW_10_2_14/src
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`
cd -
echo "python {script_loc} {thesample} -n $1 -f $2 -c $3 {gen_flag}"
time python {script_loc}  {thesample} -n $1 -f $2 -c $3 {gen_flag}
mv *.root {full_eos_out} 
'''
.format(script_loc   = script_loc, 
        thesample    = args.samples[0],
        gen_flag     = gen_flag,
        full_eos_out = full_eos_out
        )
)
subprocess.call(['chmod', '+x', bname])
    

## write the cfg for condor submission condor_multiple_readnano.cfg
with open('%s/condor_sub.cfg'%base_out, 'w') as cfg:
    cfg.write('''Universe = vanilla
Executable = {bname}
use_x509userproxy = True 
transfer_input_files = FindValueFromVectorOfBool.h
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = True
requirements = (OpSysAndVer =?= "CentOS7")
    
Log    = {base_out}/outCondor/condor_job_$(Process).log
Output = {base_out}/outCondor/condor_job_$(Process).out
Error  = {base_out}/outCondor/condor_job_$(Process).err
Arguments = $(Process) {outdir} {chan} 
Queue {njobs}
        '''.format( bname = bname, 
                     flavour = args.flavour, 
                     base_out = base_out, 
                     outdir = args.outdir, 
#                      era = args.era, 
                     chan = channel, 
                     njobs = njobs )
)    
    # submit to the queue
    print('condor_submit {base_out}/condor_sub.cfg'.format(base_out=base_out))
    if not args.test:
        os.system("condor_submit {base_out}/condor_sub.cfg".format(base_out=base_out))   


# +JobFlavour = "{flavour}"
