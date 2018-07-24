#!/usr/bin/python3

import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<annotated.genes.fasta> <clusters> <cluster_id>')
    exit()

with open(sys.argv[2]) as f:
    clusters = f.readlines()

genes = list(SeqIO.parse(sys.argv[1], 'fasta'))
cluster_of_interest = sorted(clusters[int(sys.argv[3])].strip().split())


for gene in cluster_of_interest:
    index = int(gene[5:])
    print('>{}'.format(genes[index].description))
    print('{}'.format(genes[index].seq))
