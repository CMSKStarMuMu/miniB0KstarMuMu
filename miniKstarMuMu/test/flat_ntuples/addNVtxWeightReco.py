import argparse


parser = argparse.ArgumentParser(description = "This program creates add the pileup weight to the RECO MC sample. \
The weights are taken from the corresponding GEN sample")
parser.add_argument("filenamesInGENMC",          help="Path to the input GEN MC files")
parser.add_argument("filenamesInRecoMC",         help="Path to the input Reco MC files")
parser.add_argument("-t", "--tree",              default="ntuple",            help="Name of the tree holding the variables")
# parser.add_argument("-w", "--wName",              default="weight",                 help="Name of the branch that will contain the weight")
# parser.add_argument("-c", "--cut",                default="",                       help="Cut string which is applied on number of vertices branch")
# parser.add_argument("-v", "--verbosity",          default=True,                     help="Set verbosity to [0, 1]")
args = parser.parse_args()


import os, sys, pdb
import numpy as np
import pandas, root_numpy

import root_pandas


ifile_reco = args.filenamesInRecoMC  
ifile_gen  = args.filenamesInGENMC

branch_list = ['eventN','weight']
merge_list = ['eventN']
if '2016' in ifile_gen:
    branch_list= ['eventN', 'weight', 'weightBF', 'weightGH', 'lumi']
    merge_list = ['eventN', 'lumi']

print ('loading support dataset...')
dataset_support = pandas.DataFrame(
    root_numpy.root2array(
        ifile_gen,
        'ntuple',
        branches= branch_list,
    )
)
print ('\t...done')

print ('loading dataset...')
dataset = pandas.DataFrame(
    root_numpy.root2array(
        ifile_reco,
        'ntuple',
    )
)
print ('\t...done')

dataset = pandas.merge(dataset,dataset_support[branch_list],on=merge_list, how='left')
dataset.to_root(ifile_reco.replace('.root', '_addNVtx.root'), key='ntuple', store_index=False)    