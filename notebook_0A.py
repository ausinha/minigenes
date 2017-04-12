#!/usr/bin/env python

import sys
import os
from collections import Counter
import argparse

import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import scipy.stats
import scipy.io as sio
from Bio.Seq import Seq
# from dnatools import dnatools
# %matplotlib inline

import matplotlib as mpl
mpl.use('Agg')

def reverse_complement(sequence):
    s = Seq(sequence)
    return s.reverse_complement()

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

# project_dir = os.path.dirname(args.fwd_dna)
# if not project_dir:
project_dir = '.'
if args.name:
    suffix = args.name
else:
    suffix = os.path.basename(args.fwd_rna).replace('.fq', '')

resultsdir = project_dir + '/results/' + suffix + '/'
if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)
figdir = project_dir + '/results/' + suffix + '/'
if not os.path.exists(figdir):
    os.makedirs(figdir)
datadir = project_dir + '/data/' + suffix + '/'
if not os.path.exists(datadir):
    os.makedirs(datadir)
    
# Choose if you want to actually save the plots:
SAVEFIGS = True

alt3SS_seq = '.........................GCTTGGATCTGATCTCAACAGGGT.........................'
alt3SS_tag = 'CATTACCTGC.........................'

tags = Counter()
header = {}
seq = {}
strand = {}
quality ={}
tag_seqs = {}
c = 0
p = 0
d = 0
f = {}
f[0] = open(args.fwd_dna, 'r')
f[1] = open(args.rev_dna, 'r')
while True:
    for i in range(2):
        header[i] = f[i].readline()[:-1]
        seq[i] = f[i].readline()[:-1]
        strand[i] = f[i].readline()[:-1]
        quality[i] = f[i].readline()[:-1]

    cur_tag = str(reverse_complement(seq[1]))
    if(len(header[0])==0):
        break
        
    # Check passing reads and that the sequence after the random tag matches
    # the plasmid sequence.

    # barcode is 0:20, plasmid is 20:30, (or 10:20 in rev)
    # ASK: not matching 0:10, may have duplicate barcodes
    # why not just match based on barcode
    if (cur_tag[10:20]==alt3SS_tag[:10]):
        p += 1
        #Check that the non-randomized sequences match perfectly to the reference
        if(seq[0][25:25+24]==alt3SS_seq[25:-25]):
            barcode = cur_tag[-20:]
            d+=1
            # curtag is practically barcode + fixed seq
            # dict of unique fwd read for each barcode
            try:
                tag_seqs[barcode]
            except:
                tag_seqs[barcode] = Counter()
            # only keep first 74 to match Manoj's data
            tag_seqs[barcode][seq[0][:74]]+=1

    if(c%1000000)==0:
        print c,p,d,'|',
    c+=1
print
    
for i in range(2):
    f[i].close()

# read alex data and combine with seq tags
if args.fwd_dna_my:
    if not os.path.exists(args.fwd_dna_my):
        print >>sys.stderr, 'warning: file doesnt exist', args.fwd_dna_my
    f[0] = open(args.fwd_dna_my, 'r')
    f[1] = open(args.rev_dna_my, 'r')
    c = 0
    p = 0
    d = 0
    while True:
        for i in range(2):
            header[i] = f[i].readline()[:-1]
            seq[i] = f[i].readline()[:-1]
            strand[i] = f[i].readline()[:-1]
            quality[i] = f[i].readline()[:-1]

        cur_tag = str(reverse_complement(seq[1]))
        if(len(header[0])==0):
            break
            
        # Check passing reads and that the sequence after the random tag matches
        # the plasmid sequence.

        # barcode is 0:20, plasmid is 20:30, (or 10:20 in rev)
        # ASK: not matching 0:10, may have duplicate barcodes
        # why not just match based on barcode
        # if (cur_tag[10:20]==alt3SS_tag[:10]):
        if (cur_tag[46:56]==alt3SS_tag[:10]):
            p += 1
            #Check that the non-randomized sequences match perfectly to the reference
            if(seq[0][25:25+24]==alt3SS_seq[25:-25]):
                barcode = cur_tag[-20:]
                d+=1
                # curtag is practically barcode + fixed seq
                # dict of unique fwd read for each barcode
                try:
                    tag_seqs[barcode]
                except:
                    tag_seqs[barcode] = Counter()
                # only keep first 74 to match Manoj's data
                tag_seqs[barcode][seq[0][:74]]+=1

        if(c%1000000)==0:
            print c,p,d,'|',
        c+=1
    print
        
    for i in range(2):
        f[i].close()

# for barcode in tag_seqs:
#     print barcode
#     for fwd in tag_seqs[barcode]:
#         print '\t', tag_seqs[barcode][fwd], fwd

# 755889 / 753365, less than 1% but may change the count order
# print len(tag_seqs)
# print len(set([seq[-20:] for seq in tag_seqs]))

# keep the N25A-B with highest count, incase there is spurios
ks = tag_seqs.keys()
tag_map = {}
tag_map_counts = {}
c = 0
for k in ks:
    max_seq = max(tag_seqs[k]) # Get seq
    max_seq_counts = tag_seqs[k][max_seq]
    if(max_seq_counts>=2):
        tag_map[k] = max_seq
        tag_map_counts[k] = max_seq_counts
    if(c%100000)==0:
        print c,
    c+=1
seq_series = pd.Series(tag_map)
seq_counts = pd.Series(tag_map_counts)
print

# Right now, I've kept 30nt of the barcode, even though only 20nt of that is randomized. Let's trim this to 20nt:
# TODO: no need to slice -20 coz barcode alreayd sliced
seq_series = pd.Series(dict(zip(pd.Series(seq_series.index).str.slice(-20),seq_series.values )))

filename = os.path.basename(args.fwd_dna + '.Seqs.csv')
path = datadir + '/' + filename
seq_series.name='Seq'
seq_series.index.name='Tag'
seq_series.to_csv(path, index_label='Tag', header=True)

tag2seq_dict = dict(zip(seq_series.index,np.arange(len(seq_series))))

# the first NNNs are first 25bp, then the next 25bp, then 20bp for barcode
# the 56 bp before barcode perfectly matches in rev to the DNA seq
alt3SS_full_seq = 'gtaagttatcaccttcgtggctacagagtttccttatttgtctctgttgccggcttatatggacaagcatatcacagccatttatcggagcgcctccgtacacgctattatcggacgcctcgcgagatcaatacgtataccagctgccctcgatacatgtcttggacggggtcggtgttgatatcgtatNNNNNNNNNNNNNNNNNNNNNNNNNGCTTGGATCTGATCTCAACAGGGTNNNNNNNNNNNNNNNNNNNNNNNNNatgattacacatatagacacgcgagcacccatcttttatagaatgggtagaacccgtcctaaggactcagattgagcatcgtttgcttctcgagtactacctggtacagatgtctcttcaaacaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaaggactgatagtaaggcccattacctgcNNNNNNNNNNNNNNNNNNNNGCAGAACACAGCGGTTCGACCTGCGTGATATCTCGTATGCCGTCTTCTGCTTG'
alt3SS_full_seq = alt3SS_full_seq.upper()

c = 0
header = {}
seq = {}
strand = {}
quality ={}

tag_list = []
ss_list = []

f = {}
f[0] = open(args.fwd_rna, 'r')
f[1] = open(args.rev_rna, 'r')

while True:
    for i in range(2):
        header[i] = f[i].readline()[:-1]
        seq[i] = f[i].readline()[:-1]
        strand[i] = f[i].readline()[:-1]
        quality[i] = f[i].readline()[:-1]
    if(len(header[i])==0):
        break
        #min_qual[i] = min(quality[i])
    tag = reverse_complement(seq[1][:20])

    try:
        tag_ind = tag2seq_dict[tag]
    except:
        pass
    else:
        if args.mine:
            # Check if the end of the read 100-120 matches the second exon
            # of citrine. In case of mismatches, I check for matches to 3
            # different 20nt regions.
            pos1 = 30
            pos2 = 20
            pos3 = 10
            pos4 = 0
            s_start = alt3SS_full_seq.find(seq[0][pos1:pos1+10])-pos1
            if(s_start<-pos1):
                s_start = alt3SS_full_seq.find(seq[0][pos2:pos2+10])-pos2
                if(s_start<-pos2):
                    s_start = alt3SS_full_seq.find(seq[0][pos3:pos2+10])-pos3
                    if(s_start<-pos3):
                        s_start = alt3SS_full_seq.find(seq[0][pos4:pos4+10])-pos4

                # if(s_start<-pos2):
                #     s_start = alt3SS_full_seq.find(seq[0][10:20])-10

                #     if(s_start<-10):
                #         s_start = alt3SS_full_seq.find(seq[0][0:10])-0
            if(s_start>=0):
                tag_list.append(tag_ind)
                ss_list.append(s_start)
        else:
            # Check if the end of the read 100-120 matches the second exon
            # of citrine. In case of mismatches, I check for matches to 3
            # different 20nt regions.
            s_start = alt3SS_full_seq.find(seq[0][100:120])-100
            if(s_start<-100):
                s_start = alt3SS_full_seq.find(seq[0][80:100])-80
                if(s_start<-80):
                    s_start = alt3SS_full_seq.find(seq[0][60:80])-60
            if(s_start>=0):
                tag_list.append(tag_ind)
                ss_list.append(s_start)
    if(c%1000000)==0:
        print c,
    c+=1

for i in range(2):
    f[i].close()

abc = list(np.ones_like(ss_list))+[0]
splices = {'A3SS':scipy.sparse.csr_matrix(
    (
        list(np.ones_like(ss_list))+[0],
        (
            tag_list+[len(seq_series)-1],
            ss_list+[565]
        )
    ), dtype=np.float64)}
# print
# print len(ss_list)
# print np.ones_like(ss_list)
# print list(np.ones_like(ss_list))
# print list(np.ones_like(ss_list))+[0]

# print
# print tag_list
# print tag_list+[len(seq_series)-1]
# print ss_list+[565]

'''
a[73][372] = 1
a[96][388] = 1
a[84][388] = 1
# last
a[178][565] = 0
'''

path = datadir + '/' + os.path.basename(args.fwd_rna) + '.mat'
# sio.savemat('../data/A3SS_Reads.mat',splices)
sio.savemat(path,splices)
