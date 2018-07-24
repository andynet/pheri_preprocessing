#!/usr/bin/python3

import sys
import shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv) != 5:
    print('Usage:', sys.argv[0], '<genomes.fasta> <genomes.conversion> <genes.conversion> <dir>')
    exit()

deduplicated_genomes_file = '{}/003_deduplicated.genomes.fasta'.format(sys.argv[4])
deduplicated_genomes_conversion_file = '{}/003_deduplicated.genomes.conversion'.format(sys.argv[4])
deduplicated_genes_conversion_file = '{}/003_deduplicated.genes.conversion'.format(sys.argv[4])

deduplicated_genomes = open(deduplicated_genomes_file, 'w')
shutil.copy(sys.argv[2], deduplicated_genomes_conversion_file)
shutil.copy(sys.argv[3], deduplicated_genes_conversion_file)

genomes = list(SeqIO.parse(sys.argv[1], 'fasta'))
genomes.sort(key=lambda genome: genome.seq)
genomes.sort(key=lambda genome: len(genome.seq), reverse=True)

current_genome = genomes[0]                 # type: SeqRecord
for genome in genomes:

    current_seq = current_genome.seq        # type: Seq
    current_description = current_genome.description

    if current_seq == genome.seq or current_seq.reverse_complement() == genome.seq:
        # find out lowest ID
        current_description = min(current_description, genome.description)
        # rewrite genomes conversion
        command = 'sed -i s/{}/{}/ {}'.format(genome.description, current_description, deduplicated_genomes_conversion_file)
        subprocess.call(command, shell=True)
        # rewrite genes conversion
        command = 'sed -i s/{}/{}/ {}'.format(genome.description, current_description, deduplicated_genes_conversion_file)
        subprocess.call(command, shell=True)
    else:
        deduplicated_genomes.write('>' + current_description + '\n')
        deduplicated_genomes.write(str(current_seq) + '\n')

        current_genome = genome

deduplicated_genomes.close()
