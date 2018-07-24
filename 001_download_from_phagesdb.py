#!/usr/bin/python3

from requests.exceptions import ConnectionError
from Bio import SeqIO
import requests
import time
import json
import wget
import sys
import os


def safe_get_data(feature, key, default_value='NO_DATA'):
    try:
        value = feature[key]
        if value == '':
            value = default_value
    except KeyError:
        value = default_value

    return value


def get_seq(genome):

    fasta_file = genome['fasta_file']
    fa_handle = wget.download(fasta_file)

    fa = list(SeqIO.parse(fa_handle, 'fasta'))
    if len(fa) != 1:
        print(fasta_file, 'has', len(fa), 'sequences', file=sys.stderr)

    fasta_seq = str(fa[0].seq)
    os.remove(fa_handle)

    return fasta_seq


def get_location(gene):

    start = safe_get_data(gene, 'Start')
    stop = safe_get_data(gene, 'Stop')
    orientation = safe_get_data(gene, 'Orientation')

    if orientation == 'F':
        return '[{}:{}](+)'.format(start, stop)
    if orientation == 'R':
        return '[{}:{}](-)'.format(start, stop)

    return '[{}:{}]'.format(start, stop)


def get_host(genome):

    isolation_host = safe_get_data(genome, 'isolation_host')
    if isolation_host == 'NO_DATA':
        return isolation_host
    else:
        genus = safe_get_data(isolation_host, 'genus')
        species = safe_get_data(isolation_host, 'species')
        strain = safe_get_data(isolation_host, 'strain_name')
        return ascii('{} {} {}'.format(genus, species, strain))


if len(sys.argv) != 2:
    print('Usage:', sys.argv[0], '<dir>')
    exit()

genomes_output_file = '{}/001_phagesdb.genomes.fasta'.format(sys.argv[1])
genomes_conversion_file = '{}/001_phagesdb.genomes.conversion'.format(sys.argv[1])
genomes_output = open(genomes_output_file, 'w')
genomes_conversion = open(genomes_conversion_file, 'w')

genes_output_file = '{}/001_phagesdb.genes.fasta'.format(sys.argv[1])
genes_conversion_file = '{}/001_phagesdb.genes.conversion'.format(sys.argv[1])
genes_output = open(genes_output_file, 'w')
genes_conversion = open(genes_conversion_file, 'w')

page_num = 1
genome_num = 2000000
gene_num = 20000000

while True:
    page = 'http://phagesdb.org/api/phages/?page=' + str(page_num)
    req = requests.get(page)

    if req.status_code == 404:
        print('Download finished at page', page_num-1, '.', file=sys.stderr)
        break

    if req.status_code != 200:
        time.sleep(30)
        print('\nPage', page_num, 'returned error status', req.status_code, '. Trying again.', file=sys.stderr)
        continue

    print('{} downloaded.'.format(page))
    page_num += 1

    genomes = json.loads(req.text)
    for genome in genomes['results']:

        if genome['fasta_file'] is None:
            continue

        fasta_id = 'phage{:0>7}'.format(genome_num)
        genome_num += 1

        try:
            fasta_seq = get_seq(genome)
        except UnicodeDecodeError:
            print('Cannot parse {}. Skipping...'.format(genome['fasta_file']))
            continue

        genomes_output.write('>' + fasta_id + '\n')
        genomes_output.write(fasta_seq + '\n')

        name = safe_get_data(genome, 'phage_name')
        accession = safe_get_data(genome, 'genbank_accession')
        host = get_host(genome)

        genomes_conversion.write('{}\t{}\t{}\t{}\t{}\n'.format(fasta_id, accession, name, host, 'NO_DATA'))

        try:
            page = 'http://phagesdb.org/api/genesbyphage/' + name
            req = requests.get(page)
        except ConnectionError as e:
            print('Could not get genes of {}. Skipping...'.format(name))
            print(e)
            continue

        if req.status_code != 200:
            print('\nPage', page, 'returned error status.', file=sys.stderr)
            continue

        genes = json.loads(req.text)
        for gene in genes['results']:

            gene_id = 'gene{:0>8}'.format(gene_num)
            gene_num += 1

            gene_seq = safe_get_data(gene, 'translation')

            genes_output.write('>' + gene_id + '\n')
            genes_output.write(gene_seq + '\n')

            protein_id = safe_get_data(gene, 'GeneID')
            product = safe_get_data(gene, 'Notes')
            location = get_location(gene)

            genes_conversion.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_id, fasta_id, protein_id,
                                                                     location, product, 'NO_DATA'))

genomes_output.close()
genomes_conversion.close()
genes_output.close()
genes_conversion.close()
