#!/usr/bin/python3

from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
import pandas as pd
import graphviz
import pickle
import sys
import os

if len(sys.argv) < 4:
    print('Usage:', sys.argv[0], '<spec> <matrix> <matrix>...')
    exit()

spec = sys.argv[1]
data_dir = os.path.dirname(os.path.abspath(sys.argv[2]))
file_name = os.path.basename(os.path.abspath(sys.argv[2]))
cluster_method = '.'.join(file_name.split('.')[1:-3])
mtype = file_name.split('.')[-3]

hosts = []
for arg in sys.argv[2:]:
    file_name = os.path.basename(os.path.abspath(arg))
    hosts.append(file_name.split('.')[-2])

labels = []
matrices = []

for host in hosts:
    file_name = '{}/012_matrix.{}.{}.{}.tsv'.format(data_dir, cluster_method, mtype, host)
    matrix = pd.read_csv(file_name, sep='\t', header=0, index_col=0)

    print(matrix.shape)
    if host == spec:
        labels += [1]*matrix.shape[0]
    else:
        labels += [0]*matrix.shape[0]

    matrices.append(matrix)

# merge matrices
features = pd.concat(matrices)

# build model
clf = DecisionTreeClassifier(min_impurity_split=0.03)
clf.fit(features, labels)

# visualize tree
dot_data = tree.export_graphviz(clf, out_file=None, feature_names=features.columns)
graph = graphviz.Source(dot_data)
graph.format = 'pdf'
graph.render('{}/013_tree.{}.{}.{}'.format(data_dir, cluster_method, mtype, spec))

# save model as pickle
model_file_name = '{}/013_model.{}.{}.{}.pkl'.format(data_dir, cluster_method, mtype, spec)
model_pkl = open(model_file_name, 'wb')
pickle.dump(clf, model_pkl)
model_pkl.close()
