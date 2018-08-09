#!/usr/bin/python3

import subprocess
import math
import sys
import os

if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<005_annotated.genes.conversion> <008_genes.clstr> <009_genomes.list>')
    exit()

with open(sys.argv[3]) as f:
    genome_list_lines = f.readlines()

num_of_nodes = 10
num_of_genomes = len(genome_list_lines)
genomes_per_node = math.ceil(num_of_genomes/num_of_nodes)

script_dir = sys.path[0]
data_dir = os.path.dirname(os.path.abspath(sys.argv[3]))

# <editor-fold desc=header_creating>
with open(sys.argv[2]) as f:
    cluster_num = len(f.readlines())

cluster_names = []
for i in range(cluster_num):
    cluster_names.append('Cluster_{}'.format(i))

header = '\t{}\n'.format('\t'.join(cluster_names))

with open('{}/matrix.part'.format(data_dir), 'w') as f:
    f.write(header)
# </editor-fold>

for i in range(0, num_of_genomes, genomes_per_node):
    command = '{}/009_create_matrix_from_mcl.py {} {} {} {} {} {}/matrix.part{}'.format(
        script_dir,
        sys.argv[1], sys.argv[2], sys.argv[3],
        i, min(i+genomes_per_node, num_of_genomes),
        data_dir, i//genomes_per_node)

    qsub_command = 'echo "{0}" |     ' \
                   'qsub -l thr=16  ' \
                   '     -e {1}      ' \
                   '     -o {1}      ' \
                   '     -cwd       ' \
                   '     -N matrix{2}'.format(command, '/dev/null', i//genomes_per_node)

    subprocess.call(qsub_command, shell=True)
