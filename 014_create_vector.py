#!/usr/bin/python3

import os
import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.tree import DecisionTreeClassifier

if len(sys.argv) != 7:
    print('Usage:', sys.argv[0], '<genes.blast> <cluster.tsv> <matrix> <model> <vector> <result>')
    exit()

phage_name = os.path.basename(os.path.abspath(sys.argv[1])).split('.')[0]
spec = os.path.basename(os.path.abspath(sys.argv[4])).split('.')[-2]

# create dictionary
genes_to_cluster = dict()
with open(sys.argv[2]) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    members = line.strip().split('\t')
    for member in members:
        genes_to_cluster[member] = i

# calculate vector
vector = []
with open(sys.argv[1]) as f:
    lines = f.readlines()

for line in lines:
    target_gene = line.split()[1]
    vector.append(genes_to_cluster[target_gene])

# remove duplicates and sort
vector = list(set(vector))
vector.sort()

# # write vector to file
# out = open(sys.argv[3], 'w')
# out.write('{}'.format(phage_name))
# for item in vector:
#     out.write('\t{}'.format(item))
# out.write('\n')

# load matrix and create new record for phage
matrix = pd.read_csv(sys.argv[3], sep='\t', header=0, index_col=0)
new_record = pd.DataFrame(np.zeros((1, matrix.shape[1]), dtype=int), columns=matrix.columns, index=[phage_name])

for item in vector:
    if 'Cluster_{}'.format(item) in new_record.columns:
        new_record.set_value(phage_name, 'Cluster_{}'.format(item), 1)
    else:
        print('Cluster_{} is not in the matrix.'.format(item))

new_record.to_csv(sys.argv[5], sep='\t')

# load model and classify
model_pkl = open(sys.argv[4], 'rb')
model = pickle.load(model_pkl)          # type: DecisionTreeClassifier

prediction = model.predict(new_record)

with open(sys.argv[6], 'w') as f:
    f.write('{}\t{}\t{}\n'.format(phage_name, spec, prediction))
