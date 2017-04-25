#!/usr/bin/env python
# parse positions of match from 2 bam files
# one aligned to template and the other aligned to minigenes

import argparse
from collections import Counter
import os
import sys

import pysam


def main():
    args = read_args()
    count = Counter()

    for f in args.bams:
        count = count_starts(f, count)

    for pos, num in count.most_common():
        print '{}\t{}'.format(pos, num)


def count_starts(filename, count):
    for r in pysam.Samfile(os.path.expanduser(filename), 'rb'):
        if r.is_supplementary:
            continue
        # if r.reference_start == 235:
        #     print >>sys.stderr, r.qname
        count[r.reference_start] += 1
    return count


def read_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-b', '--bams', nargs='+')
    parser.set_defaults(
    )

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
