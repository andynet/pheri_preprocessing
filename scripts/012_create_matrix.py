#!/usr/bin/python3

from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Create binary matrix')
parser.add_argument('-c', '--clusters', required=True, help='010_mcl.clusters')
parser.add_argument('-f', '--fasta', required=True, help='005_train.fna')
parser.add_argument('-o', '--output', required=True)
args = parser.parse_args()

phages = []
clusters = []

# load phage names and their number
records = list(SeqIO.parse(args.fasta, 'fasta'))
for record in records:
    phages.append(record.id)
phages.sort()

# load cluster names and their number
with open(args.clusters) as f:
    lines = f.readlines()

cluster_count = int(lines[1].split()[-1])
for i in range(cluster_count):
    clusters.append('cluster{:0>5}'.format(i+1))

# create empty SparseDataFrame
data = pd.DataFrame(np.zeros(shape=(len(phages), len(clusters))), index=phages, columns=clusters)

# for every line in clusterfile update sparse dataframe
for i in range(2, len(lines)-1):
    phage = lines[i].split()[0].split('_')[0]
    cluster = 'cluster{:0>5}'.format(int(lines[i].split()[1].strip()))

    data.at[phage, cluster] += 1

# save SparseMatrix
data.to_csv(args.output, sep='\t')
