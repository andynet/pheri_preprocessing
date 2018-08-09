#!/usr/bin/python3

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import math
import sys
import os
import re


def create_fasta_dir(directory, fasta):

    os.mkdir(directory)
    fa_list = list(SeqIO.parse(fasta, 'fasta'))

    for record in fa_list:      # type: SeqRecord

        file = '{}/{}.fa'.format(directory, record.id)
        seq_id = '>{}\n'.format(record.id)
        seq = '{0!s}\n'.format(record.seq)

        out = open(file, 'w')
        out.write(seq_id)
        out.write(seq)
        out.close()


if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<blast::tsv> <genes::fasta> <dir>')
    exit()

blast_output = sys.argv[1]
genes_fasta = sys.argv[2]

script_dir = sys.path[0]
stage_dir = sys.argv[3]
fasta_dir = '{}/fasta'.format(stage_dir)
tmp_dir = '{}/tmp'.format(stage_dir)

if not os.path.isdir(fasta_dir):
    create_fasta_dir(fasta_dir, genes_fasta)
    print('Fasta files correctly prepared.')
else:
    print('Fasta files was already prepared. Skipping...')

os.mkdir(tmp_dir)

with open(blast_output) as f:
    blast_lines = f.readlines()

# files_num is maximal number of files to be created. Bigger numbers create bigger overhead in time consumed
# by grid engine, but smaller numbers make grid engine less available for other jobs with higher priority executed
# after run of this script. For our cluster we established files_num to 80000 artificially based on number of threads
# on our cluster.
files_num = 80000
lines_num = len(blast_lines)
records_per_file = math.ceil(lines_num/files_num)
print('Final product will be {} files in format abc.'.format(math.ceil(lines_num/records_per_file)))

script_template = '#!/bin/bash\n\n'                                         \
                  '{body}\n\n'                                              \
                  'cat {list_of_needle_tsv} > {dir}/{file_num:0>5}.abc\n'   \
                  'rm {list_of_needle}; rm {list_of_needle_tsv}\n'
number = 0
commands = []
created_needle_files = []
created_needle_tsv_files = []

for i, line in enumerate(blast_lines):

    seq1 = line.split('\t')[0].strip()
    seq2 = line.split('\t')[1].strip()

    seq1_file = '{}/{}.fa'.format(fasta_dir, seq1)
    seq2_file = '{}/{}.fa'.format(fasta_dir, seq2)

    needle_out = '{}/{}_to_{}.needle'.format(tmp_dir, seq1, seq2)
    needle_tsv_out = '{}/{}_to_{}.needle.tsv'.format(tmp_dir, seq1, seq2)

    created_needle_files.append(needle_out)
    created_needle_tsv_files.append(needle_tsv_out)

    needle_command = 'needle -asequence {} -bsequence {} ' \
                     '       -outfile {}   -aformat score' \
                     '       -gapopen 10.0 -gapextend 0.5' \
                     '       -endopen 10.0 -endextend 0.5' \
                     '       -endweight Y '.format(seq1_file, seq2_file, needle_out)

    parse_command = '{}/007_parse_needle.py {}' \
                    '                       {}'.format(script_dir, needle_out, needle_tsv_out)

    commands.append(re.sub(' +', ' ', needle_command))
    commands.append(re.sub(' +', ' ', parse_command))

    if i != 0 and i % records_per_file == 0:

        script = script_template.format(body='\n'.join(commands),
                                        dir=tmp_dir, file_num=number,
                                        list_of_needle=' '.join(created_needle_files),
                                        list_of_needle_tsv=' '.join(created_needle_tsv_files))

        script_name = '{}/tmp_script_{}.sh'.format(tmp_dir, number)

        script_out = open(script_name, 'w')
        script_out.write(script)
        script_out.close()

        if number % 1000 == 0:
            qsub_command = 'echo "bash {0}; rm {0};" | '\
                           'qsub -l thr=1 '             \
                           '     -o {1} '               \
                           '     -e {1} '               \
                           '     -cwd '                 \
                           '     -N glal_{2}            '.format(script_name, '/dev/null', number)
        else:
            qsub_command = 'echo "bash {0}; rm {0};" | '\
                           'qsub -l thr=1 '             \
                           '     -o {1} '               \
                           '     -e {1} '               \
                           '     -cwd '                 \
                           '     -N glal_{2} > /dev/null'.format(script_name, '/dev/null', number)

        subprocess.call(qsub_command, shell=True)

        number += 1
        commands = []
        created_needle_files = []
        created_needle_tsv_files = []

# if i == len(blast_lines)-1
script = script_template.format(body='\n'.join(commands),
                                dir=tmp_dir, file_num=number,
                                list_of_needle=' '.join(created_needle_files),
                                list_of_needle_tsv=' '.join(created_needle_tsv_files))

script_name = '{}/tmp_script_{}.sh'.format(tmp_dir, number)

script_out = open(script_name, 'w')
script_out.write(script)
script_out.close()

qsub_command = 'echo "bash {0}; rm {0}; touch {3};" | '  \
               'qsub -l thr=1 '                         \
               '     -o {1} '                           \
               '     -e {1} '                           \
               '     -cwd '                             \
               '     -N glal_{2}            '.format(script_name, '/dev/null', number,
                                                     '{}/qsub_completed'.format(stage_dir))

subprocess.call(qsub_command, shell=True)
