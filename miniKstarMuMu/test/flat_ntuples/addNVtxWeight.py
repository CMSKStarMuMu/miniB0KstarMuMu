#!/usr/bin/env python

import argparse
import sys
import os
from array import array

"""
Setup argument parser
"""

parser = argparse.ArgumentParser(description = "This program creates a MC tree with an additional weight branch from given MC and data input files. \
The weigths are calculated from the ratio of pileup distribution in data relative to the trueNumInteraction branch MC (weight = Data/MC).")
parser.add_argument("filenamesInputData",         help="Path to the input Data files")
parser.add_argument("filenamesInputMC",           help="Path to the input MC files")
parser.add_argument("filenameOutputMC",           help="Path to the output MC file with added weight branch")
parser.add_argument("-d", "--directory",          default="tpTree",                 help="Directory in the input ROOT file which contains the input tree")
parser.add_argument("-t", "--tree",               default="fitter_tree",            help="Name of the tree holding the variables")
parser.add_argument("-b", "--branchMC",           default="trueNumInteractionsMC",  help="Branch name with number of vertices")
parser.add_argument("-B", "--branchData",         default="pileup",                 help="Branch name with number of vertices")
parser.add_argument("-w", "--wName",              default="weight",                 help="Name of the branch that will contain the weight")
parser.add_argument("-c", "--cut",                default="",                       help="Cut string which is applied on number of vertices branch")
parser.add_argument("-hf","--histogramFilename",  default="control_nVtx.png",       help="Output filename of control histogram with number of vertices")
parser.add_argument("-hr","--histogramRange",     default="100,0,100",              help="Range of control histogram given as \"bins,min,max\"")
parser.add_argument("-v", "--verbosity",          default=True,                     help="Set verbosity to [0, 1]")
args = parser.parse_args()

"""
Read input files
"""

from ROOT import * # import this here, otherwise it overwrites the argparse stuff
gROOT.SetBatch(True) # set ROOT to batch mode, this suppresses printing canvases

# Get path to tree
treePath = args.tree
# treePath = os.path.join(args.directory,args.tree)
if args.verbosity==1:
    print('Used path to tree in files:')
    print('---------------------------')
    print(treePath)
    print('')

# Setup tree chains
chainData = TFile(args.filenamesInputData,'read')
# chainData = TChain(treePath)
if args.verbosity==1:
    print('Input files data:')
    print('-----------------')
    print(args.filenamesInputData)
# for filename in args.filenamesInputData.split(' '):
#     if args.verbosity==1:
#         print(filename)
#     chainData.AddFile(filename)
if args.verbosity==1:
    print('')

chainMC = TChain(treePath)
# if args.verbosity==1:
#     print('Input files MC:')
#     print('---------------')
# for filename in args.filenamesInputMC.split(' '):
#     if args.verbosity==1:
#         print(filename)

filelistin = (os.popen('ls {FOLDER}/gen*.root'.format(FOLDER = args.filenamesInputMC)).read()).split('\n')
for filename in filelistin[:-1]:
    chainMC.AddFile(filename)

if args.verbosity==1:
    print 'n entries in gen tree: ', chainMC.GetEntries()
    print('')

"""
Make histograms of number of vertices and create control plot
"""

if args.verbosity==1:
    print('Make control histogram:')
    print('-----------------------')
c1 = TCanvas("c1", "c1")
c1.Divide(2,1)
c1.cd(1)
# chainData.Draw(args.branchData+'>>hData('+args.histogramRange+')', args.cut)
hData = chainData.Get('pileup')
chainData.Draw()
c1.cd(2)
chainMC.Draw(args.branchMC+'>>hMC('+args.histogramRange+')', args.cut)
c1.Print(args.histogramFilename)
if args.verbosity==1:
    print('')

"""
Add weight branch to MC tree chain
"""

# Calculate weights
hMC = gROOT.FindObject('hMC')
hData.Scale(1.0/hData.Integral())
hMC.Scale(1.0/hMC.Integral())

weights = [1.0]*(hData.GetNbinsX()+1)
for i in range(1, len(weights)):
    nMC = hMC.GetBinContent(i)
    nData = hData.GetBinContent(i)
    if nMC > 0:
        weights[i-1] = nData/nMC
    else:
        weights[i-1] = 1.0

# Create new file with directory
fileOutputMC = TFile.Open(args.filenameOutputMC, 'recreate')
# fileOutputMC.mkdir(args.directory).cd()
fileOutputMC.cd()

# Clone tree from input MC
progressbarWidth = 40
if args.verbosity==1:
    print('Adding weight column:')
    print('---------------------')
    sys.stdout.write('Progress: [{}]'.format('-'*progressbarWidth))
    sys.stdout.flush() # this forces to print the stdout buffer
    sys.stdout.write('\b'*(progressbarWidth+1)) # return to start of line, after '['


treeOutput = chainMC.CloneTree(0)
weight = array('f', [1.0])
treeOutput.Branch(args.wName, weight, '%s/F'%args.wName)
numEvents = chainMC.GetEntries()

# import pdb; pdb.set_trace()
for i in range(numEvents):
    chainMC.GetEntry(i)
    if int(getattr(chainMC, args.branchMC)) > len(weights):  print 'ehi, looking for ', int(getattr(chainMC, args.branchMC))
    weight[0] = weights[int(getattr(chainMC, args.branchMC))]  # NOTE we need this struct because otherwise PyROOT somehow drops objects...
    treeOutput.Fill()
#     if args.verbosity==1:
#         if i%int(numEvents/(progressbarWidth-1))==0:
#             sys.stdout.write('+')
#             sys.stdout.flush()

if args.verbosity==1:
    sys.stdout.write('\n')
    print('')

# Write file, close file and print filename
fileOutputMC.Write()
fileOutputMC.Close()
if args.verbosity==1:
    print('Output file MC:')
    print('---------------')
    print(args.filenameOutputMC)
    print('')

