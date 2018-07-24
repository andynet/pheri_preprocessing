#!/usr/bin/python3

import sys


def create_gtc_dict(file):

    with open(file) as f:
        lines = f.readlines()

    dictionary = dict()

    for line in lines:
        if line[0] == '%' or line[0] == '\n':
            continue

        phage, cluster = line.split()[0:2]
        dictionary[phage] = int(cluster)

    return dictionary


if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<009_clusters.out> <007_train.cd-hit.genes.fasta.clstr> <010_genes.tsv>')
    exit()

gene_to_cluster_dict = create_gtc_dict(sys.argv[1])

with open(sys.argv[2]) as f:
    lines = f.readlines()

cluster = None

for line in lines:
    if line[0] == '>':
        cluster = None
        continue

    phage = line.split()[2][1:-3]

    if cluster is None:
        cluster = gene_to_cluster_dict[phage]
    else:
        gene_to_cluster_dict[phage] = cluster


cluster_to_gene_dict = dict()

for key, value in gene_to_cluster_dict.items():
    if value in cluster_to_gene_dict:
        cluster_to_gene_dict[value].append(key)
    else:
        cluster_to_gene_dict[value] = [key]

with open(sys.argv[3], 'w') as f:
    for i in range(1, len(cluster_to_gene_dict)+1):
        f.write('\t'.join(cluster_to_gene_dict[i]) + '\n')

print('completed')




