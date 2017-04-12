#!/usr/bin/env python

import sys
from collections import Counter


'''search for match using rosenberg strategy, using perfect match of 10-mers'''


alt3SS_full_seq = 'gtaagttatcaccttcgtggctacagagtttccttatttgtctctgttgccggcttatatggacaagcatatcacagccatttatcggagcgcctccgtacacgctattatcggacgcctcgcgagatcaatacgtataccagctgccctcgatacatgtcttggacggggtcggtgttgatatcgtatNNNNNNNNNNNNNNNNNNNNNNNNNGCTTGGATCTGATCTCAACAGGGTNNNNNNNNNNNNNNNNNNNNNNNNNatgattacacatatagacacgcgagcacccatcttttatagaatgggtagaacccgtcctaaggactcagattgagcatcgtttgcttctcgagtactacctggtacagatgtctcttcaaacaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaaggactgatagtaaggcccattacctgcNNNNNNNNNNNNNNNNNNNNGCAGAACACAGCGGTTCGACCTGCGTGATATCTCGTATGCCGTCTTCTGCTTG'
alt3SS_full_seq = alt3SS_full_seq.upper()

count = Counter()

while True:
    header = sys.stdin.readline().strip()
    seq = sys.stdin.readline().strip()
    strand = sys.stdin.readline().strip()
    quality = sys.stdin.readline().strip()

    if not header:
        break

    # Check if the end of the read 100-120 matches the second exon
    # of citrine. In case of mismatches, I check for matches to 3
    # different 20nt regions.
    pos1 = 30
    pos2 = 20
    pos3 = 10
    pos4 = 0
    s_start = alt3SS_full_seq.find(seq[pos1:pos1+10])-pos1
    if(s_start<-pos1):
        s_start = alt3SS_full_seq.find(seq[pos2:pos2+10])-pos2
        if(s_start<-pos2):
            s_start = alt3SS_full_seq.find(seq[pos3:pos2+10])-pos3
            if(s_start<-pos3):
                s_start = alt3SS_full_seq.find(seq[pos4:pos4+10])-pos4

        # if(s_start<-pos2):
        #     s_start = alt3SS_full_seq.find(seq[0][10:20])-10

        #     if(s_start<-10):
        #         s_start = alt3SS_full_seq.find(seq[0][0:10])-0
    if(s_start>=0):
        count[s_start] += 1

for pos, num in count.most_common():
    print '{}\t{}'.format(pos, num)
