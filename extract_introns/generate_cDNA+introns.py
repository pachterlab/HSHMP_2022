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
        fh.write(data)
        # fh.write('\n'.join([data[i:i+80] for i in range(0, len(data), 80)]))
        fh.write('\n')

def collapse_data(intervals, sequence, strand):
    data = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
    data = collapse_N(data.upper())
    if strand == '-':
        data = reverse_complement(data)
    return data

def process_gene(gene, scaffolds, out, nascent, mature):
    try:
        sequence = scaffolds[gene['scaffold']]
    except KeyError as _:
        print(f'Scaffold {gene["scaffold"]} not found.')
        return

    print(f'processing {gene["name"]}')

    if not nascent and not mature:
        nascent = True
        mature = True

    if nascent:
        # Find start of first exon and end of last exon in gene
        end = 0
        start = float('inf')
        for tr, data in gene['trs'].items():
            start = min(start, min(e['start'] for e in data['exons']))
            end = max(end, max(e['end'] for e in data['exons']))
        seq = sequence[start:end]
        header = f'>{gene["name"]} nascent_transcript chromosome:GRCh38:{gene["scaffold"]}:{gene["start"]}:{gene["end"]} gene:{gene["name"]} strand{gene["strand"]}\n'
        out.write(header)
        # out.write('\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]))
        out.write(seq)
        out.write('\n')

    for tr, data in gene['trs'].items():
        exons = data['exons']
        # TODO:
        # Decide how to handle single-exon isoforms

        # Nascent transcript
        # seq = sequence[gene['start']:gene['end']]
        # # if gene['strand'] == '-':
            # # seq = reverse_complement(seq)

        # if nascent:
            # header = f'>{tr}.N nascent_transcript chromosome:GRCh38:{gene["scaffold"]}:{gene["start"]}:{gene["end"]} gene:{gene["name"]} strand{gene["strand"]}\n'
            # out.write(header)
            # # out.write('\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]))
            # out.write(seq)
            # out.write('\n')

        # Mature transcript
        if mature:
            ivs = sorted([(e['start'], e['end']) for e in exons], key=lambda e: e[0])
            ivs = list(merge_intervals(ivs))
            header = f'>{tr} mature_transcript chromosome:GRCh38:{gene["scaffold"]}:{ivs[0][0]}:{ivs[-1][1]} gene:{gene["name"]} strand{gene["strand"]}\n'
            seq = collapse_data(ivs, sequence, '+')
            out.write(header)
            out.write('\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]))
            out.write('\n')

def parse_gtf(path, scaffolds, out, nascent, mature):
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

    fh = open(out, 'w')
    for _, gene in genes.items():
        process_gene(gene, scaffolds, fh, nascent, mature)
    fh.close()

def generate_cDNA_introns(gtf_path, fasta_path, out='.', nascent=False, mature=False):
    scaffolds = parse_fasta(fasta_path)
    parse_gtf(gtf_path, scaffolds, out, nascent, mature)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nascent', action='store_true')
    parser.add_argument('--mature', action='store_true')
    parser.add_argument('--gtf', type=str, help='Path to GTF file')
    parser.add_argument('--fa', type=str, help='Path to fasta file')
    parser.add_argument('--out', type=str, help='Path to output file')
    args = parser.parse_args()

    generate_cDNA_introns(args.gtf, args.fa, args.out, args.nascent, args.mature)
