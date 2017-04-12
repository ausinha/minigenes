#!/usr/bin/env python

import argparse

import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import scipy.stats
import os
import scipy.io as sio

import matplotlib as mpl
mpl.use('Agg')

import dnatools
from plot_tools import simpleaxis, plot_splicing_histogram
fsize=14
try:
    from pylab import *
    # Plotting Params:
    rc('mathtext', default='regular')
    fsize=14
except ImportError, e:
    print >>sys.stderr, e

parser = argparse.ArgumentParser(description='')
parser.add_argument('--name')
parser.add_argument('--fwd-dna')
parser.add_argument('--rev-dna')
parser.add_argument('--fwd-dna-my')
parser.add_argument('--rev-dna-my')
parser.add_argument('--fwd-rna')
parser.add_argument('--rev-rna')
parser.add_argument('--mine', action='store_true')

parser.set_defaults(
)
args = parser.parse_args()

project_dir = '.'
if args.name:
    suffix = args.name
else:
    suffix = os.path.basename(args.fwd_rna).replace('.fq', '')

resultsdir = project_dir + '/results/' + suffix + '/'
if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)
# figdir = '../figures/N1_Library_Statistics/'
figdir = project_dir + '/figures/' + suffix + '/'
if not os.path.exists(figdir):
    os.makedirs(figdir)
datadir = project_dir + '/data/' + suffix + '/'
if not os.path.exists(datadir):
    os.makedirs(datadir)
    
#Choose if you want to actually save the plots:
SAVEFIGS = True

path = datadir + os.path.basename(args.fwd_rna) + '.mat'
print path
data = sio.loadmat(path)
# A5SS_data = data['A5SS']
# print data

# A3SS
A3SS_data = data['A3SS']
# Only look at SA_1 usage:
A3SS_data = np.array(A3SS_data[:,235].todense()).reshape(-1)/np.array(A3SS_data.sum(axis=1),dtype=np.float64).reshape(-1)
# Get minigenes with reads
A3SS_nn = find(pd.notnull(A3SS_data))
A3SS_data = A3SS_data[A3SS_nn]
filename = os.path.basename(args.fwd_dna + '.Seqs.csv')
path = datadir + filename
print path
A3SS_seqs = pd.read_csv(path,index_col=0).Seq[A3SS_nn]

# Make a web logo for splice acceptors in the A3SS library

#Get each randomized region
r1 = pd.Series(A3SS_seqs).str.slice(0,25)
r2 = pd.Series(A3SS_seqs).str.slice(-25)
Y = data['A3SS'][A3SS_nn]

base_freq_all = {}
for b in dnatools.bases:
    base_freq_all[b] = r1.str.count(b).sum() + r2.str.count(b).sum()
base_freq_all = pd.Series(base_freq_all)

#Find all new SA in the two randomized regions . Then calculate the frequency of each base at each position
sa_base_freqs = {}
for p in range(15,26):
    inds = find(np.array(Y[:,189+24+25+p].todense())>0)
    for b in range(0,min(p+5,25)):
        try:
            sa_base_freqs[b-p] = sa_base_freqs[b-p].add(r2[inds].groupby(r2[inds].str.slice(b,b+1)).size(),fill_value=0)
        except:
            sa_base_freqs[b-p] = r2[inds].groupby(r2[inds].str.slice(b,b+1)).size()
for p in range(15,26):
    inds = find(np.array(Y[:,189+p].todense())>0)
    for b in range(0,min(p+5,25)):
        try:
            sa_base_freqs[b-p] = sa_base_freqs[b-p].add(r2[inds].groupby(r1[inds].str.slice(b,b+1)).size(),fill_value=0)
        except:
            sa_base_freqs[b-p] = r2[inds].groupby(r1[inds].str.slice(b,b+1)).size()
sa_base_freqs = pd.DataFrame(sa_base_freqs)

sa_to_random_ratio = (sa_base_freqs.fillna(0).T/base_freq_all).T
normalized_base_ratio = (sa_to_random_ratio/sa_to_random_ratio.T.sum(axis=1))

base_prob_dict = normalized_base_ratio.to_dict()
sampled_bases = {}
for pos in range(-25,5):
    sampled_bases[pos] = []
    for b in 'ACGT':
        sampled_bases[pos] += [b]*int(round(base_prob_dict[pos][b]*1000000))
    sampled_bases[pos] = np.array(sampled_bases[pos])
def make_sa():
    sa = ''
    for i in range(-25,5):
        sa += np.random.choice(sampled_bases[i])
    return sa

f = open(resultsdir+'Generated_SA.txt','w')
for i in range(10000):
    f.write(make_sa()+'\n')
f.close()

# weblogo=matplotlib.image.imread(figdir+'Library_SA_Normalized_Frequencies.png')
# import ipdb; ipdb.set_trace()
# print weblogo
# print type(weblogo)

