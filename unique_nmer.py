#!/usr/bin/env python

import argparse


def main():
    args = read_args()
    seq = open('../docs/minigene.txt').read().strip()

    k = args.kmer
    for i in range(0, len(seq) - k + 1):
        sub = seq[i:i+k]
        pos = seq.find(sub, i + 1)
        if pos >= 0 and sub != 'N' * k:
            print i, pos, seq[i-2:i], seq[pos-2:pos], sub


def read_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', '--param')
    parser.add_argument('-k', '--kmer', type=int)
    parser.set_defaults(
        kmer=15,
    )

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
