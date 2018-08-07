#!/usr/bin/python3

from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
import pandas as pd
import argparse
import graphviz
import pickle
import os


def load_lines(filename):
    with open(filename) as f:
        lines = f.readlines()
    return lines


def create_phage_to_genus_dict(lines):
    result = dict()
    for line in lines:
        phage_id = line.split()[0]
        genera = line.split()[1].strip()
        if phage_id not in result:
            result[phage_id] = set()
        result[phage_id].add(genera)

    return result


parser = argparse.ArgumentParser(description='Create decision tree models and their pdf representations')
parser.add_argument('-m', '--matrix', required=True)
parser.add_argument('-g', '--genera', required=True)
parser.add_argument('-i', '--genera_of_interest', required=True)
parser.add_argument('-d', '--data_dir', required=True)
args = parser.parse_args()

features = pd.read_csv(args.matrix, sep='\t', header=0, index_col=0)
phage_to_genus = create_phage_to_genus_dict(load_lines(args.genera))
genera_of_interest = [item.strip() for item in load_lines(args.genera_of_interest)]
phages = features.index.values

data_dir = args.data_dir
file_name = os.path.basename(os.path.abspath(args.matrix))

cluster_method = file_name.split('.')[1]
matrix_type = file_name.split('.')[2]

for genus_of_interest in genera_of_interest:
    labels = []

    for phage in phages:
        if genus_of_interest in phage_to_genus[phage]:
            labels.append(1)
        else:
            labels.append(0)

    # build model
    clf = DecisionTreeClassifier(max_leaf_nodes=20)
    clf.fit(features, labels)

    # visualize tree
    dot_data = tree.export_graphviz(clf, out_file=None, feature_names=features.columns)
    graph = graphviz.Source(dot_data)
    graph.format = 'pdf'
    graph.render('{}/013_{}.tree.{}.{}'.format(data_dir, genus_of_interest, cluster_method, matrix_type))

    # save model as pickle
    model_file_name = '{}/013_{}.model.{}.{}.pkl'.format(data_dir, genus_of_interest, cluster_method, matrix_type)
    model_pkl = open(model_file_name, 'wb')
    pickle.dump(clf, model_pkl)
    model_pkl.close()