#!/usr/bin/python3

import os
import sys
import pandas as pd

if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<011_matrix.{cluster_method}.{type}.tsv> <004_hosts> <HOSTS>')
    exit()

matrix_file = sys.argv[1]
matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)

# hosts loading
hosts_dict = dict()

with open(sys.argv[2]) as f:
    for line in f:
        phage = line.split()[0]
        string = line.split()[1]

        hosts_dict[phage] = string

# calculating how to split matrix
queries = sys.argv[3].split(',')
query = dict()

for q in queries:
    query[q] = []

for i in range(matrix.shape[0]):
    name = matrix.iloc[i].name

    qnum = 0
    for q in queries:
        if q in hosts_dict[name]:
            if qnum == 0:
                query[q].append(name)
            else:
                print('Warning: {} {} could be in more groups.'.format(name, hosts_dict[name]))
            qnum += 1

if sum([len(query[q]) for q in queries]) < matrix.shape[0]:
    print('Some phages may not have host in queries')

# creating split matrices
dir_name = os.path.dirname(matrix_file)
methods = '.'.join(matrix_file.split('.')[1:-1])

for q in queries:

    new_matrix = matrix.loc[query[q]]
    matrix_out = '{}/012_matrix.{}.{}.tsv'.format(dir_name, methods, q)
    new_matrix.to_csv(matrix_out, sep='\t')

    hosts_out = '{}/012_hosts.{}.{}'.format(dir_name, methods, q)
    hosts_out_content = ''
    for record in query[q]:
        hosts_out_content += '{}\t{}\n'.format(record, hosts_dict[record])
    with open(hosts_out, 'w') as f:
        f.write(hosts_out_content)
