#!/usr/bin/env python

import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# the first NNNs are first 25bp, then the next 25bp, then 20bp for barcode
# the 56 bp before barcode perfectly matches in rev to the DNA seq
alt3SS_full_seq = 'gtaagttatcaccttcgtggctacagagtttccttatttgtctctgttgccggcttatatggacaagcatatcacagccatttatcggagcgcctccgtacacgctattatcggacgcctcgcgagatcaatacgtataccagctgccctcgatacatgtcttggacggggtcggtgttgatatcgtatNNNNNNNNNNNNNNNNNNNNNNNNNGCTTGGATCTGATCTCAACAGGGTNNNNNNNNNNNNNNNNNNNNNNNNNatgattacacatatagacacgcgagcacccatcttttatagaatgggtagaacccgtcctaaggactcagattgagcatcgtttgcttctcgagtactacctggtacagatgtctcttcaaacaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaaggactgatagtaaggcccattacctgcNNNNNNNNNNNNNNNNNNNNGCAGAACACAGCGGTTCGACCTGCGTGATATCTCGTATGCCGTCTTCTGCTTG'

part1 = alt3SS_full_seq[0:189]
part2 = alt3SS_full_seq[214:238]
part3 = alt3SS_full_seq[263:612]
part4 = alt3SS_full_seq[632:]

in_file = 'data/MPR34AB_001.10bp/A3SS_dna_R1.fq.Seqs.csv'
minigenes = []
with open(in_file) as fh:
    fh.readline()
    for k, line in enumerate(fh):
        if k % 10000 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()

        barcode, rest = line.strip().split(',')
        n25a = rest[:25]
        n25b = rest[-25:]

        minigene = alt3SS_full_seq[:]
        minigene = ''.join([part1, n25a, part2, n25b, part3, barcode, part4])
        seq = Seq(minigene.lower())

        record = SeqRecord(seq,'chr' + str(k + 1), '', barcode)
        minigenes.append(record)

        # print '>chr' + str(k + 1)
        # cSeqIO.write(sequences, fasta, 'fasta')
        # print record.format('fasta')
        # if k > 10:
            # break

# with open('minigene.fa', 'w') as fasta:
SeqIO.write(minigenes, 'minigene.fa', 'fasta')
