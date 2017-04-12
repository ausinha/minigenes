#!/usr/bin/env python
# trim variable 6-11 bp based on motif

import sys
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--fwd')
parser.add_argument('--rev')
args = parser.parse_args()

ifwd = open(args.fwd)
irev = open(args.rev)

ofwd = open(args.fwd.replace('.fq', '.trim.fq'), 'w')
orev = open(args.rev.replace('.fq', '.trim.fq'), 'w')

while True:
    header = ifwd.readline().strip()
    if not header:
        break
    seq = ifwd.readline().strip()
    misc = ifwd.readline().strip()
    quality = ifwd.readline().strip()
    
    header_rev = irev.readline().strip()
    seq_rev = irev.readline().strip()
    misc_rev = irev.readline().strip()
    quality_rev = irev.readline().strip()
    
    # print both fwd and rev if only motif found, keeps pair in sync
    pos = seq.find('GATCCGCCACAACATCGAG')
    if pos >= 6 and pos <= 11:
        print >>ofwd, header
        print >>ofwd, seq[pos+19:]
        print >>ofwd, misc
        print >>ofwd, quality[pos+19:]

        print >>orev, header_rev
        print >>orev, seq_rev
        print >>orev, misc_rev
        print >>orev, quality_rev
