#!/usr/bin/env python3

"""
Author:
    Kristján Eldjárn Hjörleifsson
    keldjarn@caltech.edu

Usage:
./count_kmers.py --gtf annotation.gtf3 --fa scaffolds.fasta

"""

from operator import itemgetter
from itertools import chain
import argparse
import gzip

from utils import (reverse_complement, parse_rest, parse_fasta, collapse_N,
                   merge_intervals, interval_diff)

def collapse_data(intervals, sequence, strand):
    data = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
    data = collapse_N(data.upper())
    if strand == '-':
        data = reverse_complement(data)
    return data

def process_gene(gene, scaffolds):
    try:
        sequence = scaffolds[gene['scaffold']]
    except KeyError as _:
        print(f'Scaffold {gene["scaffold"]} not found.')
        return

    print(f'processing {gene["name"]}')

    nascent = set()
    mature = set()
    ambiguous = set()
    for tr, data in gene['trs'].items():
        exons = data['exons']

        ivs = sorted([(e['start'], e['end']) for e in exons], key=lambda e: e[0])
        ivs = list(merge_intervals(ivs))
        if len(ivs) > 1:

            # Get k-mers containing exon-exon boundaries
            for idx in range(len(ivs) - 1):
                seq = sequence[max(ivs[idx][0], ivs[idx][1]-30) : ivs[idx][1]] + sequence[ivs[idx+1][0] : min(ivs[idx+1][1], ivs[idx+1][0]+30)]
                if len(seq) > 30:
                    for start in range(len(seq) - 30):
                        mature.add(seq[start:start+31])

        # Get all interior k-mers of exons
        for iv in ivs:
            seq = sequence[iv[0]:iv[1]]
            if len(seq) > 30:
                for start in range(len(seq) - 30):
                    ambiguous.add(seq[start:start+31])

        # All other k-mers are nascent
        seq = sequence[gene['start']:gene['end']]
        if len(seq) > 30:
            for start in range(len(seq) - 30):
                kmer = seq[start:start+31]
                if kmer in mature:
                    ambiguous.add(kmer)
                elif kmer not in ambiguous:
                    nascent.add(kmer)

    return nascent, mature, ambiguous


def parse_gtf(path, scaffolds):
    # Who the heck came up with this hecking file format?
    genes = {}
    # Apparently, gtf files are not necessarily ordered by gene, so we cannot
    # do this in a single pass-through
    with gzip.open(path, 'r') as fh:
        for l in fh:

            line = l.decode('utf-8')
            # Skip comments
            if line.startswith('#'):
                continue

            data = line.split('\t')
            fields = {
                'scaffold': data[0],
                'feature': data[2],
                'start': int(data[3]) - 1, # Scaffold is 0-indexed
                'end': int(data[4]) - 1,
                'strand': data[6],
                'rest': parse_rest(data[8])
            }
            gene_id = fields['rest']['gene_id']

            if gene_id not in genes:
                genes[gene_id] = {
                    'name': gene_id,
                    'trs': {},
                    'scaffold': fields['scaffold'],
                    'strand': fields['strand']
                }

            if fields['feature'] in ['exon',
                                     'UTR',
                                     'start_codon',
                                     'stop_codon',
                                     'five_prime_utr',
                                     'three_prime_utr',
                                     'CDS']:
                tr = fields['rest']['transcript_id']
                if tr not in genes[gene_id]['trs']:
                    genes[gene_id]['trs'][tr] = {
                        'name': tr,
                        'exons': [fields]
                    }
                else:
                    genes[gene_id]['trs'][tr]['exons'].append(fields)
            elif fields['feature'] == 'gene':
                genes[gene_id]['start'] = fields['start']
                genes[gene_id]['end'] = fields['end']

    nascent = set()
    mature = set()
    ambiguous = set()
    for _, gene in genes.items():
        n, m, a = process_gene(gene, scaffolds)
        nascent = {*nascent, *n}
        mature = {*mature, *m}
        ambiguous = {*ambiguous, *a}

    nascent_count = 0
    mature_count = 0
    ambiguous_count = 0
    for m in mature:

        revcomp = reverse_complement(m)
        if m in ambiguous or revcomp in ambiguous:
            continue

        if m in nascent or revcomp in nascent:
            ambiguous.add(m)
            continue

        mature_count += 1

    for n in nascent:

        revcomp = reverse_complement(n)
        if n in ambiguous or revcomp in ambiguous:
            continue

        if n in mature or revcomp in mature:
            ambiguous.add(n)
            continue

        nascent_count += 1

    print(f'Nascent: {nascent_count}, mature: {mature_count}, ambiguous: {len(ambiguous)}')

def count_kmers(gtf_path, fasta_path):
    scaffolds = parse_fasta(fasta_path)
    parse_gtf(gtf_path, scaffolds)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', type=str, help='Path to GTF file')
    parser.add_argument('--fa', type=str, help='Path to fasta file')
    args = parser.parse_args()

    count_kmers(args.gtf, args.fa)
