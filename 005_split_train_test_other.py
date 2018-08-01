#!/usr/bin/python3

from random import sample
from Bio import SeqIO
from math import floor
import argparse
import os

parser = argparse.ArgumentParser(description='Split phages to train and test datasets.')
parser.add_argument('-f', '--fasta', required=True, help='Fasta file containing genomes.')
parser.add_argument('-g', '--genera', required=True, help='Data about phage hosts')
parser.add_argument('-l', '--list', required=True, help='File containing genera of interests - one per line.')
parser.add_argument('-r', '--ratio', type=float, default=0.8)
args = parser.parse_args()

with open(args.list) as f:
    genera = f.readlines()

with open(args.genera) as f:
    lines = f.readlines()

train_set = set()
test_set = set()

for category in genera:
    category = category.strip()

    tmp_set = set()
    for line in lines:
        phage = line.split()[0].strip()
        genus = line.split()[1].strip()

        if genus == category:
            tmp_set.add(phage)

    tmp_train_set = set(sample(tmp_set, floor(args.ratio * len(tmp_set))))
    tmp_test_set = tmp_set - tmp_train_set

    train_set = train_set | tmp_train_set       # union
    test_set = test_set | tmp_test_set          # union

print('Train set: {}\nTest set: {}\n'.format(len(train_set), len(test_set)))

intersection = train_set & test_set
if len(intersection) != 0:
    print('Sets are not disjoint sets. {} is in both sets.'.format(intersection))

data_dir = os.path.dirname(os.path.abspath(args.fasta))

# split genome.fasta
genome_fa = list(SeqIO.parse(args.fasta, 'fasta'))

train_out = open('{}/005_train.fna'.format(data_dir), 'w')
test_out = open('{}/005_test.fna'.format(data_dir), 'w')
other_out = open('{}/005_other.fna'.format(data_dir), 'w')

for record in genome_fa:        # type: Bio.Seq

    if record.id in train_set:
        train_out.write('>{}\n{}\n'.format(record.id, record.seq))
    elif record.id in test_set:
        test_out.write('>{}\n{}\n'.format(record.id, record.seq))
    else:
        other_out.write('>{}\n{}\n'.format(record.id, record.seq))

train_out.close()
test_out.close()
other_out.close()
