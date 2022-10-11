#!/usr/bin/env python3

"""
Author:
    Kristján Eldjárn Hjörleifsson
    keldjarn@caltech.edu

Usage:
./extract_introns.py --gtf annotation.gtf3 --fa scaffolds.fasta -out output_directory [--union] [--diff]

"""

from operator import itemgetter
from itertools import chain
import argparse
import gzip

from utils import (reverse_complement, parse_rest, parse_fasta, collapse_N,
                   merge_intervals, interval_diff)

def write_gene(gene, data, out):
    if len(data) == 0:
        return
    with open(path, 'a') as fh:
        fh.write(f'>{tr}\n')
        fh.write('\n'.join([data[i:i+80] for i in range(0, len(data), 80)]))
        fh.write('\n')

def collapse_data(intervals, sequence, strand):
    data = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
    data = collapse_N(data.upper())
    if strand == '-':
        data = reverse_complement(data)
    return data

def process_gene(gene, scaffolds, out, union, diff):
    try:
        sequence = scaffolds[gene['scaffold']]
    except KeyError as _:
        print(f'Scaffold {gene["scaffold"]} not found.')
        return
    transcripts = {
        tr: [(e['start'], e['end']) for e in exons['exons']]\
             for tr, exons in gene['trs'].items()
    }
    ivs = list(chain(*transcripts.values()))
    start = min(ivs, key=lambda i: i[0])[0]
    end = max(ivs, key=lambda i: i[1])[1]

    exon_union = list(merge_intervals(ivs))
    itrs = [(e1[1], e2[0]) for _, exons in transcripts.items()\
                           for e1, e2 in zip(exons[:-1], exons[1:])]

    # Look for retained introns at beginning/end of gene
    for _, exons in transcripts.items():
        if exons[0][0] > start:
            itrs.append((start, exons[0][0]))
        if exons[-1][1] < end:
            itrs.append((exons[-1][1], end))
    itrs = sorted(itrs, key=itemgetter(0))

    if union:
        itrs = list(merge_intervals(itrs))

        if diff:
            itrs = interval_diff(itrs, exon_union)

    for i, iv in enumerate(itrs):
        out.write(f'>{gene["name"]}.{"union." if union else ""}{"exondiff." if diff else""}tr.{i}\n')
        seq = sequence[iv[0]:iv[1]]
        if gene['strand'] == '-':
            seq = reverse_complement(seq)
        out.write('\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]))
        out.write('\n')

def parse_gtf(path, scaffolds, out, union, diff):
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

    fh = open(out, 'w')
    for _, gene in genes.items():
        process_gene(gene, scaffolds, fh, union, diff)
    fh.close()

def extract_introns(gtf_path, fasta_path, out='.', union=False, diff=False):
    scaffolds = parse_fasta(fasta_path)
    parse_gtf(gtf_path, scaffolds, out, union, diff)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--union', action='store_true')
    parser.add_argument('--diff', action='store_true')
    parser.add_argument('--gtf', type=str, help='Path to GTF file')
    parser.add_argument('--fa', type=str, help='Path to fasta file')
    parser.add_argument('--out', type=str, help='Path to output file')
    args = parser.parse_args()

    extract_introns(args.gtf, args.fa, args.out, args.union, args.diff)
