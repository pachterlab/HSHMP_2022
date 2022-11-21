#!/usr/bin/env python3

from collections.abc import Iterator
import sys
from bisect import bisect_left
import argparse

def revcomp(string):
    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'M': 'K',
        'K': 'M',
        'Y': 'R',
        'R': 'Y',
        'S': 'S',
        'W': 'W',
        'N': 'N'
    }
    return ''.join(list(map(lambda c: comp[c], string))[::-1])

def sequence_set_diff(set_a, set_b):

    set_b.sort()

    for seq in set_a:

        idx = bisect_left(set_b, seq)

        if idx >= len(set_b) or set_b[idx] != seq:

            yield seq.rstrip()

def parse_mature_reads(path):

    capture = False

    with open(path, 'r') as fh:

        for count, line in enumerate(fh):

            if count % 4 == 0 and 'mature' in line:

                capture = True

            elif count % 4 == 1 and capture:

                capture = False
                yield(line.rstrip())

def count_ambiguous_reads(path: str, mature: set) -> int:

    capture = False
    ambiguous = 0

    with open(path, 'r') as fh:

        for count, line in enumerate(fh):

            if count % 4 == 0 and 'nascent' in line:

                capture = True

            elif count % 4 == 1 and capture:

                capture = False
                if line.rstrip() in mature:

                    ambiguous += 1

    return ambiguous

def parse_txome(path: str) -> str:
    # Generate a searchable string from transcriptome
    txome_db = ''

    with open(path) as fh:

        for line in fh:

            if line[0] == '>':
                # Need to demarcate the transcripts so as not to align
                # across the boundaries

                txome_db += '|'

            else:
                txome_db += line.rstrip('\n')

    return txome_db

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--mapped', type=str, help='Mapped reads')
    parser.add_argument('--mapped-diff', type=str, help='Reads to be subtracted from --mapped')
    parser.add_argument('--reads', type=str, help='Simulated reads file')
    parser.add_argument('--txome', type=str, help='The reference transcriptome. Warning: compute-heavy')
    args = parser.parse_args()

    txome_db = ''
    if args.txome is not None:
        # Generate a searchable string from transcriptome

        print(f'Parsing transcriptome from {args.txome}')
        txome_db = parse_txome(args.txome)

    if args.mapped_diff is not None:
        # Only process reads that are in --mapped and not in --mapped-diff

        print(f"Subtracting {args.mapped_diff} from {args.mapped}.")
        diff = sequence_set_diff(open(args.mapped).readlines(),
                                 open(args.mapped_diff).readlines())
        print(f"Finding mature reads in diff")
    else:
        # Process all reads in --mapped

        diff = (line.rstrip() for line in open(args.mapped))
        print("Finding mature reads in mapped reads")

    # Number of mature reads in diff
    mature_count = 0
    diff_size = 0
    # Mature reads ground truth
    mature = list(parse_mature_reads(args.reads))
    mature_size = len(mature)
    mature = set(mature)
    # Should be empty if off-list works (fingers crossed)
    ambiguous_nascent_mapped = []

    checked = 0
    for read in diff:

        diff_size += 1
        # Look for read as-is
        if read in mature:
            mature_count += 1

        elif args.txome is not None: # and diff_size % 100 == 0:

            checked += 1
            if read in txome_db or revcomp(read) in txome_db:

                ambiguous_nascent_mapped.append(read)

            else:
                print(read)


    ambiguous_in_ground_truth = count_ambiguous_reads(args.reads, mature)

    print(f"Out of {diff_size} reads in mapped reads, {mature_count} of which were mature")
    print(f"Simulated reads contain {mature_size} mature reads and {ambiguous_in_ground_truth} nascent reads are duplicated in the mature reads")
    print(f"Checked {checked} mapped nascent transcripts and {len(ambiguous_nascent_mapped)} were ambiguous")
    # print(f'Out of the {diff_size - mature_count} nascent reads that got mapped, {len(ambiguous_nascent_mapped)} were ambiguous')
