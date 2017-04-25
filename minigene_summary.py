#!/usr/bin/env python

import argparse
import os

import pysam


def main():
    args = read_args()
    minigenes = {}

    for bam in args.bams:
        for r in pysam.Samfile(os.path.expanduser(bam), 'rb'):
            if r.is_supplementary:
                continue
            splice_pos = r.reference_start

            if r.reference_name not in minigenes:
                minigenes[r.reference_name] = {
                    'transcripts': [],
                    'splice_positions': {},
                }
            if splice_pos not in \
                    minigenes[r.reference_name]['splice_positions']:
                minigenes[r.reference_name]['splice_positions'][splice_pos] \
                    = 0

            minigenes[r.reference_name]['transcripts'].append(r.qname)
            minigenes[r.reference_name]['splice_positions'][splice_pos] += 1

    for ref, info in minigenes.iteritems():
        positions = ['{}:{}'.format(pos, count)
                     for pos, count in info['splice_positions'].iteritems()]
        print '\t'.join([
            ref,
            ';'.join(info['transcripts']),
            ';'.join(positions),
        ])


def read_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-b', '--bams', nargs='+')
    parser.set_defaults(
    )

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
