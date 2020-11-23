import networkx as nx
import pickle
import random
import rdflib
import re
import pandas as pd
from rwr import partition

pkl_triple_path = '/Users/michael/PycharmProjects/VariantRankingPheKnowLator/Data/PheKnowLator_Subclass_RelationsOnly_NotClosed_NoOWL_Triples_Identifiers.txt'

# remove snps
with open('Edgelists/pkl_triples_snp_free.tsv', 'w') as file:
    for line in open(pkl_triple_path, 'r'):
        if not re.match('.*\/snp\/.*', line) and not re.match('.*http://purl.obolibrary.org/obo/BFO_0000001.*',
                                                              line) and not re.match(
            '.*http://purl.obolibrary.org/obo/doid.*', line):
            file.write(line)

# convert gene symbols
id_symbol_map = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in
                 open('/Users/michael/PycharmProjects/VariantRankingPheKnowLator/Data/gene_id2symbol.txt', 'r')}

f = 0
nf = 0
nf_set = set()
# rename nodes to use gene symbol rather than id
with open('Edgelists/pkl_edgelist_gene_symbols_no_snps.tsv', 'w') as out:
    for line in open('Edgelists/pkl_triples_snp_free.tsv'):
        line = line.strip()
        edge = line.split('\t')
        outline = ''
        # convert first node
        if re.match('.*gene.*', edge[0]):
            gene = re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[0])
            if gene in id_symbol_map:

                outline = id_symbol_map[gene] + '\t' + edge[1] + '\t'
            else:
                outline = edge[0] + '\t' + edge[1] + '\t'
        else:
            outline = edge[0] + '\t' + edge[1] + '\t'
        # convert second node
        if re.match('.*gene.*', edge[2]):
            gene = re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[2])
            if gene in id_symbol_map:
                outline = id_symbol_map[gene] + '\n'
            else:
                outline += edge[2] + '\n'
        else:
            outline += edge[2] + '\n'
        out.write(outline)

G = nx.Graph()
# for each disease gene set
for line in open('Edgelists/pkl_edgelist_gene_symbols_no_snps.tsv'):
    edge = line.strip().split('\t')
    if len(edge) != 3:
        continue
    G.add_edge(edge[0], edge[2], predicate=edge[1])

gene_sets = 'DisGeNET_genesets.txt'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]

disease_genes = disease_gene_sets[list(disease_gene_sets.keys())[0]]
nodes = partition(disease_genes, 2)
seeds, targets = nodes[0], nodes[1]


counts = {}
for disease in disease_gene_sets.keys():
    disease_genes = disease_gene_sets[disease]
    nodes = partition(disease_genes, 2)
    seeds, targets = nodes[0], nodes[1]
    for i,seed in enumerate(seeds):
        print(seed,str(i),':',str(len(seeds)))
        for target in targets:
            try:
                path = nx.shortest_path(G, seed, target)
            except nx.exception.NodeNotFound:
                print('Not in graph',seed,'or',target)
            for i in range(len(path) - 1):
                p = G[path[i]][path[i + 1]]['predicate']
                if p in counts:
                    counts[p] += 1
                else:
                    counts[p] = 1
    # break