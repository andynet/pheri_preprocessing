#!/usr/bin/python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Split big fasta file to records.')
parser.add_argument('-i', '--fasta', required=True)
parser.add_argument('-o', '--output_dir', required=True)
args = parser.parse_args()

fasta = list(SeqIO.parse(args.fasta, 'fasta'))

for record in fasta:    # type: SeqRecord

    with open('{}/{}'.format(args.output_dir, record.id), 'w') as f:
        f.write('>{}\n{}\n'.format(record.id, record.seq))
