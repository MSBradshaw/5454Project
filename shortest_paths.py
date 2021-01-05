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
d_keys = list(disease_gene_sets.keys())
for i,disease in enumerate(disease_gene_sets.keys()):
    print(i)
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
                'http://purl.obolibrary.org/obo/CLO_0000167': 'has disease',
                'http://purl.obolibrary.org/obo/RO_0000086': 'has quality'}

named_counts = {name_mapping[x]: counts[x] for x in counts.keys()}

import matplotlib.pyplot as plt

df = pd.DataFrame({'name': list(named_counts.keys()), 'count': list(named_counts.values())})
df = df.sort_values('name')
# plt.hist(df['name'], df['count'])

import seaborn as sns

g = sns.barplot(df['name'], df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('predicate_counts.png')
plt.show()


"""
Randomly shuffle the seeds and targets
"""

from random import randrange
counts = {}
d_keys = list(disease_gene_sets.keys())

for i, disease in enumerate(disease_gene_sets.keys()):
    print(i)
    shuffled_index = None
    # keep trying to get a shuffled index until it is not the current index
    while shuffled_index is None or shuffled_index == i:
        shuffled_index = randrange(len(d_keys))

    disease_genes = disease_gene_sets[disease]
    nodes = partition(disease_genes, 2)
    seeds, targets = nodes[0], nodes[1]

    shuffled_disease_genes = disease_gene_sets[d_keys[shuffled_index]]
    shuffled_nodes = partition(shuffled_disease_genes, 2)
    shuffled_seeds, shuffled_targets = shuffled_nodes[0], shuffled_nodes[1]

    # use the shuffled seeds but the original targets
    for i, seed in enumerate(shuffled_seeds):
        # print(seed, str(i), ':', str(len(seeds)))
        for target in targets:
            try:
                path = nx.shortest_path(G, seed, target)
            except nx.exception.NodeNotFound:
                pass
                # print('Not in graph', seed, 'or', target)
            for i in range(len(path) - 1):
                p = G[path[i]][path[i + 1]]['predicate']
                if p in counts:
                    counts[p] += 1
                else:
                    counts[p] = 1
    # break
pickle.dump(counts, open('shuffled_counts.pickle', 'wb'))
counts = pickle.load(open('shuffled_counts.pickle', 'rb'))

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
                'http://purl.obolibrary.org/obo/CLO_0000167': 'has disease',
                'http://purl.obolibrary.org/obo/RO_0000086': 'has quality',
                'http://purl.obolibrary.org/obo/so#part_of': 'part of',
                'http://purl.obolibrary.org/obo/BFO_0000051': 'has part'}

named_counts = {name_mapping[x]: counts[x] for x in counts.keys()}

import matplotlib.pyplot as plt

df = pd.DataFrame({'name': list(named_counts.keys()), 'count': list(named_counts.values())})
df = df.sort_values('name')
# plt.hist(df['name'], df['count'])

import seaborn as sns

g = sns.barplot(df['name'], df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('shuffled_predicate_counts.png')
plt.show()

shuffled_counts = pickle.load(open('shuffled_counts.pickle', 'rb'))
counts = pickle.load(open('counts.pickle', 'rb'))


named_counts = {name_mapping[x]: counts[x] for x in counts.keys()}
shuffled_named_counts = {name_mapping[x]: shuffled_counts[x] for x in shuffled_counts.keys()}

diff_counts = {}
scaled_diff_counts = {}
for key in set(list(named_counts.keys()) +list(shuffled_named_counts.keys())):
    if key in named_counts:
        c = named_counts[key]
    else:
        c = 0
    if key in shuffled_named_counts:
        sc = shuffled_named_counts[key]
    else:
        sc = 0
    diff_counts[key] = sc - c
    if c is not 0:
        scaled_diff_counts[key] = (sc - c) / c
    elif sc < 0:
        scaled_diff_counts[key] = -1
    else:
        scaled_diff_counts[key] = 1

df = pd.DataFrame({'Predicate': list(scaled_diff_counts.keys()), 'Proportional Change': list(scaled_diff_counts.values())})
df = df.sort_values('Predicate')
df.plot(kind='bar',y='Proportional Change',x='Predicate')
plt.savefig('proportional_change_shuffled_predicate_counts.png')
plt.tight_layout()
plt.show()

df['Proportional Change'] = [abs(x) for x in df['Proportional Change']]
df.to_csv('proportional_change_weights.tsv',header=False,index=False,sep='\t')