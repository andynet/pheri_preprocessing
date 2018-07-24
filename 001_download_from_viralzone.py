#!/usr/bin/python3

from http.client import IncompleteRead
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError
from Bio import Entrez
from Bio import SeqIO
import time
import sys


def safe_get_qualifier(feature, key, default_value='NO_DATA'):
    try:
        return ';'.join(feature.qualifiers[key])
    except KeyError:
        return default_value


if len(sys.argv) != 2:
    print('Usage:', sys.argv[0], '<dir>')
    exit()

genomes_output_file = '{}/001_viralzone.genomes.fasta'.format(sys.argv[1])
genomes_conversion_file = '{}/001_viralzone.genomes.conversion'.format(sys.argv[1])
genomes_output = open(genomes_output_file, 'w')
genomes_conversion = open(genomes_conversion_file, 'w')

genes_output_file = '{}/001_viralzone.genes.fasta'.format(sys.argv[1])
genes_conversion_file = '{}/001_viralzone.genes.conversion'.format(sys.argv[1])
genes_output = open(genes_output_file, 'w')
genes_conversion = open(genes_conversion_file, 'w')

# query = '"Viruses"[Organism] AND srcdb_refseq[PROP] NOT wgs[PROP] NOT "cellular organisms"[Organism] ' \
#         'NOT AC_000001[PACC] : AC_999999[PACC] AND "vhost bacteria"[Filter] AND (viruses[filter] ' \
#         'AND biomol_genomic[PROP] AND ("10000"[SLEN] : "100000000"[SLEN]))'

query = '''
        (phage[Title] OR bacteriophage[Title]) 
        AND (complete sequence[Title] OR complete genome[Title]) 
        AND (viruses[filter] AND biomol_genomic[PROP] AND refseq[filter])
        '''

Entrez.email = 'andrejbalaz001@gmail.com'

handle = Entrez.esearch(db="nucleotide", term=query, idtype='acc', retmax=100000, usehistory='y')
record = Entrez.read(handle)
handle.close()

webenv = record["WebEnv"]
query_key = record["QueryKey"]
batch = 100
genome_num = 1000000
gene_num = 10000000

i = 0
while i < len(record['IdList']):

    print('Downloaded', i, 'sequences.')
    sys.stdout.flush()

    fa_handle = None
    gb_handle = None

    try:
        fa_handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", idtype="acc",
                                  retstart=i, retmax=batch, webenv=webenv, query_key=query_key)
        gb_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                  retstart=i, retmax=batch, webenv=webenv, query_key=query_key)
    except HTTPError as e:
        print(e)
        time.sleep(10)
        if fa_handle is not None:
            fa_handle.close()
        if gb_handle is not None:
            gb_handle.close()
        continue

    try:
        fa_list = list(SeqIO.parse(fa_handle, 'fasta'))
        gb_list = list(SeqIO.parse(gb_handle, 'genbank'))
    except IncompleteRead as e:
        print(e)
        time.sleep(10)
        continue
    except ValueError as e:
        print(e)
        time.sleep(10)
        continue
    finally:
        fa_handle.close()
        gb_handle.close()

    for j in range(0, len(fa_list)):

        try:
            fa = fa_list[j]                 # type: SeqRecord
            gb = gb_list[j]                 # type: SeqRecord
        except IndexError as e:
            print(e)
            i -= batch
            break

        fasta_id = 'phage{:0>7}'.format(genome_num)
        genome_num += 1

        fasta_seq = str(fa.seq)

        genomes_output.write('>' + fasta_id + '\n')
        genomes_output.write(fasta_seq + '\n')

        for feature in gb.features:
            if feature.type == 'source':

                organism = safe_get_qualifier(feature, 'organism')
                host = safe_get_qualifier(feature, 'host')
                lab_host = safe_get_qualifier(feature, 'lab_host')

                genomes_conversion.write('{}\t{}\t{}\t{}\t{}\n'.format(fasta_id, gb.id, organism, host, lab_host))

            if feature.type == 'CDS':

                gene_id = 'gene{:0>8}'.format(gene_num)
                gene_num += 1

                gene_seq = safe_get_qualifier(feature, 'translation')

                genes_output.write('>' + gene_id + '\n')
                genes_output.write(gene_seq + '\n')

                protein_id = safe_get_qualifier(feature, 'protein_id')
                product = safe_get_qualifier(feature, 'product')
                note = safe_get_qualifier(feature, 'note')

                genes_conversion.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_id, fasta_id, protein_id,
                                                                         feature.location, product, note))

    i += batch

genomes_output.close()
genomes_conversion.close()
genes_output.close()
genes_conversion.close()
