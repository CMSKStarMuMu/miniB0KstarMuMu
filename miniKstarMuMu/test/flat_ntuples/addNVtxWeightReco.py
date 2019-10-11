import argparse


parser = argparse.ArgumentParser(description = "This program creates add the pileup weight to the RECO MC sample. \
The weights are taken from the corresponding GEN sample")
parser.add_argument("filenamesInGENMC",          help="Path to the input GEN MC files")
parser.add_argument("filenamesInRecoMC",         help="Path to the input Reco MC files")
# parser.add_argument("filenamesOutRecoMC",        help="Path to the output Reco MC file with added weight branch")
# parser.add_argument("-d", "--directory",          default="tpTree",                 help="Directory in the input ROOT file which contains the input tree")
parser.add_argument("-t", "--tree",               default="fitter_tree",            help="Name of the tree holding the variables")
# parser.add_argument("-w", "--wName",              default="weight",                 help="Name of the branch that will contain the weight")
# parser.add_argument("-c", "--cut",                default="",                       help="Cut string which is applied on number of vertices branch")
# parser.add_argument("-v", "--verbosity",          default=True,                     help="Set verbosity to [0, 1]")
args = parser.parse_args()


import os, sys
import numpy as np
import pandas, root_numpy

import root_pandas


ifile_reco = args.filenamesInRecoMC  
ifile_gen  = args.filenamesInGENMC

print 'loading support dataset...'
dataset_support = pandas.DataFrame(
    root_numpy.root2array(
        ifile_gen,
        'ntuple',
        branches= ['eventN','weight'],##feat_names + additional,
    )
)
print '\t...done'


print 'loading dataset...'
dataset = pandas.DataFrame(
    root_numpy.root2array(
        ifile_reco,
        'ntuple',
#           branches =  additional + ['l1_12_5','l1_11_4']
    )
)
print '\t...done'

dataset_support['weight'] = dataset_support.weight

# dataset['weight'] = dataset_support.loc[dataset_support.eventN.isin(dataset['eventN'])].eventN
# dataset['weight'] = dataset.loc[dataset.eventN.isin(dataset_support['eventN'])].eventN

dataset.merge(dataset_support)
dataset = pandas.merge(dataset,dataset_support[['eventN','weight']],on='eventN', how='left')
# import pdb; pdb.set_trace()

dataset.to_root(ifile_reco.replace('.root', '_addNVtx.root'), key='ntuple', store_index=False)

# samples = [
#            'data_LMNR',
#            'data_Charmonium',
#            'MC_LMNR', 
#            'MC_JPSI', 
# #            'MC_BuJpsiK', 
# #            'MC_LambdaB', 
#           ]
# 
# for str_file in samples:
#     for i in range(0,11):  
#       
#         tag = '__preliminary_noiso_secondOptimization_conda_'+ str(i)
#     #     tag = '_for_subsample_newnewSB_yesMassFlatten_yesSubSample_yesKstarMass_' + str(i)
#          
#         ifile = 'sub_samples/sample_%s_events_addL1_%s.root'%(str_file, str(i))  
#       
#       
#         ## add branch for BDT 
#         print 'computing probabilities...'
#         bdt_prob_array = classifier.predict_proba(dataset_support[feat_names])[:,1]
#         print '\t...done'
#         
#         print 'adding new column to the dataset...'
#         dataset['bdt_prob'] = bdt_prob_array
#         print '\t...done'	
#         
# #         import pdb; pdb.set_trace()
#         # https://github.com/scikit-hep/root_pandas
#         dataset.to_root(ifile.replace('.root', '_addBDT.root'), key='ntuple', store_index=False)
# #         root_numpy.array2root(dataset, 'selected_tree.root', 'ntuple')
#         # ifile.GetName().replace('.root', '_%s.root'%tag)
    