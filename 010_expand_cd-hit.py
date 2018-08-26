#!/usr/bin/python3

import argparse


parser = argparse.ArgumentParser(description='Reverse deduplicating with cd-hit')
parser.add_argument('clusters')
parser.add_argument('cdhit_clusters')
parser.add_argument('out_file')
args = parser.parse_args()

with open(args.clusters) as f:
    lines = f.readlines()

gene_to_cluster = dict()
for i in range(2, len(lines)-1):
    gene = lines[i].split()[0]
    cluster = int(lines[i].split()[1].strip())

    gene_to_cluster[gene] = cluster

lines = lines[0:-1]

with open(args.cdhit_clusters) as f:
    cdhit_clstr = f.readlines()

cluster = None
for line in cdhit_clstr:
    if line[0] == '>':
        cluster = None
    else:
        gene = line.split()[2][1:-3]
        if gene in gene_to_cluster:
            cluster = gene_to_cluster[gene]
        else:
            lines.append('{}\t{}\n'.format(gene, cluster))

lines.append('\n')
with open(args.out_file, 'w') as f:
    f.writelines(lines)
