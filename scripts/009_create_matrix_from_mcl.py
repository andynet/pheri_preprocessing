#!/usr/bin/python3

import multiprocessing
import sys


def get_genes_to_clusters(lines):

    result = dict()
    for i, line in enumerate(lines):
        records = line.split()
        for record in records:
            result[record] = 'cluster{}'.format(i)

    return result


def get_phage_genes(phage, genes):

    result = []
    for gene in genes:
        gene_id = gene.split('\t')[0]
        phage_id = gene.split('\t')[1]

        if phage_id == phage:
            result.append(gene_id)

    return result


def get_vector(phage_number):

    global genes_conv
    global genes_to_clusters
    global genomes_list
    global number_of_clusters

    phage = genomes_list[phage_number].strip()
    genes = get_phage_genes(phage, genes_conv)

    active_clusters = set()
    for gene in genes:
        try:
            active_clusters.add(genes_to_clusters[gene])
        except KeyError:
            print('Gene', gene, 'has no cluster.')

    vector = []
    for i in range(0, number_of_clusters):
        if 'cluster{}'.format(i) in active_clusters:
            vector.append('1')
        else:
            vector.append('0')

    return '{}\t{}\n'.format(phage, '\t'.join(vector))


def init(genes_conv_pointer, genes_to_cluster_pointer, phage_list_pointer, number_of_clusters_pointer):

    global genes_conv
    global genes_to_clusters
    global genomes_list
    global number_of_clusters

    genes_conv = genes_conv_pointer
    genes_to_clusters = genes_to_cluster_pointer
    genomes_list = phage_list_pointer
    number_of_clusters = number_of_clusters_pointer

# <editor-fold desc="start processes">
if len(sys.argv) != 7:
    print('Usage:', sys.argv[0], '<005_annotated.genes.conversion> <008_genes.clstr> <009_genomes.list> ',
                                 '<start> <stop> <out>')
    exit()

with open(sys.argv[1]) as f:
    genes_conv = f.readlines()

with open(sys.argv[2]) as f:
    genes_clstr = f.readlines()
    number_of_clusters = len(genes_clstr)
    genes_to_clusters = get_genes_to_clusters(genes_clstr)

with open(sys.argv[3]) as f:
    genomes_list = f.readlines()

start = int(sys.argv[4])
stop = int(sys.argv[5])

pool = multiprocessing.Pool(16, initializer=init, initargs=(genes_conv, genes_to_clusters,
                                                            genomes_list, number_of_clusters))
result_array = pool.map(get_vector, range(start, stop))

with open(sys.argv[6], 'w') as f:
    f.write(''.join(result_array))
# </editor-fold>
