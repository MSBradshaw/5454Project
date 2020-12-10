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
    for i, seed in enumerate(seeds):
        print(seed, str(i), ':', str(len(seeds)))
        for target in targets:
            try:
                path = nx.shortest_path(G, seed, target)
            except nx.exception.NodeNotFound:
                print('Not in graph', seed, 'or', target)
            for i in range(len(path) - 1):
                p = G[path[i]][path[i + 1]]['predicate']
                if p in counts:
                    counts[p] += 1
                else:
                    counts[p] = 1
    # break
pickle.dump(counts, open('counts.pickle', 'wb'))
counts = pickle.load(open('counts.pickle', 'rb'))

name_mapping = {'http://purl.obolibrary.org/obo/RO_0000056': 'participates_in',
                'http://purl.obolibrary.org/obo/RO_0002160': 'evolutionarily related to',
                'http://www.w3.org/2000/01/rdf-schema#subClassOf': 'sub class of',
                'http://purl.obolibrary.org/obo/RO_0003302': 'causes or contributes to condition',
                'http://purl.obolibrary.org/obo/RO_0002511': 'Inverse of transcribed from',
                'http://purl.obolibrary.org/obo/RO_0001025': 'located in',
                'http://purl.obolibrary.org/obo/RO_0002205': 'has gene product',
                'http://purl.obolibrary.org/obo/RO_0002606': 'is substance that treats',
                'http://purl.obolibrary.org/obo/RO_0002434': 'interacts with',
                'http://purl.obolibrary.org/obo/RO_0002200': 'has phenotype',
                'http://purl.obolibrary.org/obo/BFO_0000050': 'part of',
                'http://purl.obolibrary.org/obo/RO_0000052': 'inheres in',
                'http://purl.obolibrary.org/obo/CLO_0000179': ' is disease model for',
                'http://purl.obolibrary.org/obo/RO_0002435': 'genetically interacts with',
                'http://purl.obolibrary.org/obo/RO_0004019': 'disease has basis in',
                'http://purl.obolibrary.org/obo/CLO_0000015': 'derives from patient having disease',
                'http://purl.obolibrary.org/obo/RO_0002436': 'molecularly interacts with',
                'http://purl.obolibrary.org/obo/RO_0002573': 'has modifier',
                'http://purl.obolibrary.org/obo/IDO_0000664': 'has_material_basis_in',
                'http://purl.obolibrary.org/obo/so#has_part': 'has part',
                'http://purl.obolibrary.org/obo/so#non_functional_homolog_of': 'non_functional_homolog_of',
                'http://purl.obolibrary.org/obo/so#member_of': 'member of',
                'http://purl.obolibrary.org/obo/RO_0002314': 'inheres in part of',
                'http://purl.obolibrary.org/obo/CLO_0000167': 'has disease'}

named_counts = {name_mapping[x]: counts[x] for x in counts.keys()}

import matplotlib.pyplot as plt

df = {'name': list(named_counts.keys()), 'count': list(named_counts.values())}

plt.hist(df['name'], df['count'])

import seaborn as sns

g = sns.barplot(df['name'], df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('predicate_counts.png')
plt.show()


