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
import matplotlib.pyplot as plt

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
parser.add_argument('-d', '--data')

parser.set_defaults(
)
args = parser.parse_args()

data = sio.loadmat(args.data)
A3SS_data = data['A3SS']
# A5SS_seqs = pd.read_csv('../data/A5SS_Seqs.csv',index_col=0)


# Size of Libraries:
# A5SS_reads_per_plasmid = np.array(A5SS_data.sum(1)).flatten()
A3SS_reads_per_plasmid = np.array(A3SS_data.sum(1)).flatten()
plt.hist(A3SS_reads_per_plasmid, 50, normed=1, log=True)
plt.title(os.path.basename(args.data).replace('.fq.mat', ''))
plt.xlabel('Number of reads')
plt.ylabel('Frequency (Log)')
plt.savefig('reads_per_plasmid.' + os.path.basename(args.data) + '.png')
# import ipdb; ipdb.set_trace()

# A5SS_nonzero_reads = A5SS_reads_per_plasmid>0
A3SS_nonzero_reads = A3SS_reads_per_plasmid>0

# A5SS_num_plasmids = shape(A5SS_data)[0]
A3SS_num_plasmids = shape(A3SS_data)[0]

# print 'A5SS:', sum(A5SS_nonzero_reads),'/',A5SS_num_plasmids,sum(A5SS_nonzero_reads)/float(A5SS_num_plasmids)
print 'A3SS:', sum(A3SS_nonzero_reads),'/',A3SS_num_plasmids,sum(A3SS_nonzero_reads)/float(A3SS_num_plasmids)

print 'Total Reads - A3SS: %i' %A3SS_data.sum()

# Average reads per plasmid
A3SS_total_reads = float(A3SS_data.sum())
# A5SS_total_reads = float(A5SS_data.sum())
# print A5SS_total_reads/A5SS_num_plasmids
print 'Avg read per plasmid', A3SS_total_reads/A3SS_num_plasmids

nn = A3SS_reads_per_plasmid>0
# import ipdb; ipdb.set_trace()
sys.exit(0)

for k in range(684):
    try:
        num_reads = int(A3SS_data[:,k].sum())
        if num_reads:
            print '{}\t{}'.format(k, num_reads)
    except IndexError:
        continue

sys.exit(0)
# import ipdb; ipdb.set_trace()

# Plotting Histograms of Splice Site Usage in Each Library
A3_isoforms = {}
A3_isoforms['SA_1']=np.array(A3SS_data[:,235].todense()).reshape(-1).astype(np.float64)[nn]/A3SS_reads_per_plasmid[nn]
A3_isoforms['SA_2']=np.array(A3SS_data[:,388].todense()).reshape(-1).astype(np.float64)[nn]/A3SS_reads_per_plasmid[nn]
A3_isoforms['SA_CRYPT']=np.array(A3SS_data[:,388-16].todense()).reshape(-1).astype(np.float64)[nn]/A3SS_reads_per_plasmid[nn]
A3_isoforms['SA_NEW']=np.array(A3SS_data[:,388-19-25:388-19].sum(axis=1)+A3SS_data[:,388+3:388+28].sum(axis=1)).reshape(-1).astype(np.float64)[nn]/A3SS_reads_per_plasmid[nn]
A3_isoforms['No_SA']=np.array(A3SS_data[:,0].todense()).reshape(-1).astype(np.float64)[nn]/A3SS_reads_per_plasmid[nn]

nn = A3SS_reads_per_plasmid>0
plot_splicing_histogram('SA_1','$SA_1$',A3_isoforms['SA_1'],
                        y_max=.050,
                        y_step=.010,
                        smallplot=True,
                        ylabel=True,
                        scale='M',
                        figdir=figdir,
                        save_plot=SAVEFIGS)
plot_splicing_histogram('SA_2','$SA_2$',A3_isoforms['SA_2'],
                        y_max=.100,
                        y_step=.02,
                        smallplot=True,
                        ylabel=False,
                        mean_pos=0.35,
                        smallplotpos=0.33,
                        scale='M',
                        figdir=figdir,
                        save_plot=SAVEFIGS)
plot_splicing_histogram('SA_NEW','$SA_{NEW}$',A3_isoforms['SA_NEW'],
                        y_max=.005,
                        y_step=.001,
                        smallplot=True,
                        ylabel=False,
                        scale='M',
                        figdir=figdir,
                        save_plot=SAVEFIGS)
plot_splicing_histogram('SA_CRYPT','$SA_{CRYPT}$',A3_isoforms['SA_CRYPT'],
                        y_max=.050,
                        y_step=.010,
                        smallplot=True,
                        ylabel=False,
                        scale='M',
                        figdir=figdir,
                        save_plot=SAVEFIGS)
plot_splicing_histogram('No_SA','No SA',A3_isoforms['No_SA'],
                        y_max=.050,
                        y_step=.010,
                        smallplot=True,
                        ylabel=False,
                        scale='M',
                        figdir=figdir,
                        save_plot=SAVEFIGS)

# A3SS:
A3_fraction_by_SA = np.array(A3SS_data.sum(axis=0),dtype=np.float).flatten()/A3SS_data.sum()
cryptic_SAs = (A3_fraction_by_SA>=1e-7) & (A3_fraction_by_SA<=5e-3)
print 'Number of different unique cryptic SA positions:',sum(cryptic_SAs)
print 'Percentage of reads mapping to these cryptic SAs:',A3_fraction_by_SA[cryptic_SAs].sum()

multiple_read_minigenes = find(A3SS_reads_per_plasmid>1)
np.array(A3SS_data[multiple_read_minigenes,0].todense()).flatten()

multiple_read_minigenes = find(A3SS_reads_per_plasmid>1)
SA1_multiple_reads = np.array(A3SS_data[multiple_read_minigenes,235].todense(),
                              dtype=np.float).flatten()/\
                     A3SS_reads_per_plasmid[multiple_read_minigenes]
    
print float(sum(SA1_multiple_reads==1))/len(SA1_multiple_reads),\
      sum(SA1_multiple_reads==1),'/',len(SA1_multiple_reads)
