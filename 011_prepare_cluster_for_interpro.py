#!/usr/bin/python3


import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Get genes from cluster')
parser.add_argument('fasta')
parser.add_argument('clusters')
parser.add_argument('cluster_id', type=int)
args = parser.parse_args()


with open(args.clusters) as f:
    lines = f.readlines()

gene_to_cluster = dict()
for i in range(2, len(lines)-1):
    gene = lines[i].split()[0]
    cluster = int(lines[i].split()[1].strip())

    gene_to_cluster[gene] = cluster

genes = list(SeqIO.parse(args.fasta, 'fasta'))
for gene in genes:
    if gene_to_cluster[gene.description] == args.cluster_id:
        print('>{}'.format(gene.description))
        print('{}'.format(gene.seq))
