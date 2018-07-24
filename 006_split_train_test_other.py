#!/usr/bin/python3

import os
import sys
import random
from Bio import SeqIO

if len(sys.argv) != 8:
    print('Usage:', sys.argv[0], '005_annotated.genes.conversion 005_annotated.genes.fasta',
                                 '003_deduplicated.genomes.conversion 003_deduplicated.genomes.fasta',
                                 '101_hosts',
                                 'mycobac,strepto,escheri,gordoni,pseudom,arthrob,lactoco,staphyl',
                                 '0.8')
    exit()

train_set = set()
other_set = set()

with open(sys.argv[5]) as f:
    lines = f.readlines()

selected_hosts = sys.argv[6].split(',')
train_ratio = float(sys.argv[7])

for line in lines:
    phage, host_string = line.split()

    select = False
    for host in selected_hosts:
        if host in host_string:
            select = True

    if select:
        train_set.add(phage)
    else:
        other_set.add(phage)

train_set_count = int(len(train_set) * train_ratio)

train_list = list(train_set)
random.shuffle(train_list)

train_set = set(train_list[:train_set_count+1])
test_set = set(train_list[train_set_count+1:])

print('Train set: {}, Test set: {}, Other set: {}'.format(len(train_set), len(test_set), len(other_set)))

data_dir = os.path.dirname(os.path.abspath(sys.argv[1]))

# split genome.conversion
with open(sys.argv[3]) as f:
    genome_conv = f.readlines()

suffix = '.genomes.conversion'
train_out = open('{}/006_train{}'.format(data_dir, suffix), 'w')
test_out = open('{}/006_test{}'.format(data_dir, suffix), 'w')
other_out = open('{}/006_other{}'.format(data_dir, suffix), 'w')

for record in genome_conv:
    phage = record.split('\t')[0]

    if phage in train_set:
        train_out.write(record)
    if phage in test_set:
        test_out.write(record)
    if phage in other_set:
        other_out.write(record)

train_out.close()
test_out.close()
other_out.close()

# split genome.fasta
genome_fa = list(SeqIO.parse(sys.argv[4], 'fasta'))

suffix = '.genomes.fasta'
train_out = open('{}/006_train{}'.format(data_dir, suffix), 'w')
test_out = open('{}/006_test{}'.format(data_dir, suffix), 'w')
other_out = open('{}/006_other{}'.format(data_dir, suffix), 'w')

for record in genome_fa:        # type: Bio.Seq

    if record.id in train_set:
        train_out.write('>{}\n{}\n'.format(record.id, record.seq))
    if record.id in test_set:
        test_out.write('>{}\n{}\n'.format(record.id, record.seq))
    if record.id in other_set:
        other_out.write('>{}\n{}\n'.format(record.id, record.seq))

train_out.close()
test_out.close()
other_out.close()

# split genes.conversion
with open(sys.argv[1]) as f:
    genes_conv = f.readlines()

suffix = '.genes.conversion'
train_out = open('{}/006_train{}'.format(data_dir, suffix), 'w')
test_out = open('{}/006_test{}'.format(data_dir, suffix), 'w')
other_out = open('{}/006_other{}'.format(data_dir, suffix), 'w')

train_set_genes, test_set_genes, other_set_genes = set(), set(), set()

for record in genes_conv:
    gene, phage = record.split('\t')[0:2]

    if phage in train_set:
        train_out.write(record)
        train_set_genes.add(gene)
    if phage in test_set:
        test_out.write(record)
        test_set_genes.add(gene)
    if phage in other_set:
        other_out.write(record)
        other_set_genes.add(gene)

train_out.close()
test_out.close()
other_out.close()

# split genes.fasta
genes_fa = list(SeqIO.parse(sys.argv[2], 'fasta'))

suffix = '.genes.fasta'
train_out = open('{}/006_train{}'.format(data_dir, suffix), 'w')
test_out = open('{}/006_test{}'.format(data_dir, suffix), 'w')
other_out = open('{}/006_other{}'.format(data_dir, suffix), 'w')

for record in genes_fa:

    if record.id in train_set_genes:
        train_out.write('>{}\n{}\n'.format(record.id, record.seq))
    if record.id in test_set_genes:
        test_out.write('>{}\n{}\n'.format(record.id, record.seq))
    if record.id in other_set_genes:
        other_out.write('>{}\n{}\n'.format(record.id, record.seq))

train_out.close()
test_out.close()
other_out.close()
