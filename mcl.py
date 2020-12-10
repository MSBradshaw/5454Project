import networkx as nx
import markov_clustering as mc
from rwr import create_numbered_edgelist
from diamond import partition
import pandas as pd
import sys


#
# gene_sets = 'DisGeNET_genesets.txt'
# results_name = 'test'
# results_df = None
# scores = {'500 count': [], '500 %': [], '100 count': [], '100 %': [], '50 count': [], '50 %': [], '25 count': [],
#               '25 %': [], '10 count': [], '10 %': [], 'disease': []}
#
# disease_gene_sets = {}
# for line in open(gene_sets, 'r'):
#     row = line.strip().split('\t')
#     # the first item is the disease common name, everything else is genes
#     disease_gene_sets[row[0]] = row[1:]
# for disease in disease_gene_sets.keys()
#     disease_genes = disease_gene_sets[disease]
#
#     edgelist = 'Edgelists/pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv'
#     num_edgelist_name = edgelist.split('.')[0] + '_numbered' + edgelist.split('.')[1]
#     mapping = create_numbered_edgelist(edgelist,num_edgelist_name)
#
#     inverse_mapping = { mapping[key]:key for key in mapping.keys()}
#
#     G = nx.read_edgelist(num_edgelist_name)
#
#     matrix = nx.to_scipy_sparse_matrix(G)
#
#     result = mc.run_mcl(matrix)           # run MCL with default parameters
#     clusters = mc.get_clusters(result)    # get clusters
#     for i in range(50):
#         # split the disease set randomly in 2
#         nodes = partition(disease_genes, 2)
#         seeds, targets = nodes[0], nodes[1]
#
#         # rank genes for based on number of targets in it's cluster
#         nodes = []
#         target_count = []
#
#         # convert seeds to their numerical equivalent or -1 is they are not in the network
#         numbered_seeds = [mapping[x] if x in mapping else -1 for x in seeds if x in mapping]
#         numbered_targets = [mapping[x] if x in mapping else -1 for x in targets if x in mapping]
#
#         for c in clusters:
#             num_seeds_in_cluster = sum(t in c for t in numbered_seeds)
#             for n in c:
#                 # make sure that we only include nodes that are genes
#                 if 'http' not in inverse_mapping[n]:
#                     nodes.append(n)
#                     target_count.append(num_seeds_in_cluster)
#
#         df = pd.DataFrame({'nodes':nodes,'target_count':target_count})
#         df = df.sort_values('target_count',ascending=False)
#
#         num_targets_in_network = sum(str(mapping[t]) in list(G.nodes) for t in targets if t in mapping)
#         top_xs = [500, 100, 50, 25, 10, 0]
#         for i in range(len(top_xs[:-1])):
#             print(set(df.iloc[top_xs[i + 1]:top_xs[i], 0]))
#             top_overlap = len(list(set(df.iloc[top_xs[i + 1]:top_xs[i], 0]) & set(numbered_targets)))
#             scores[str(top_xs[i]) + ' count'].append((top_overlap))
#             scores[str(top_xs[i]) + ' %'].append(top_overlap / num_targets_in_network)
#         scores['disease'].append(disease)
#         single_df = pd.DataFrame(scores)
#         # get the mean of each of the iterations for this disease
#         means = single_df.mean(axis=0)
#         print(means)
#         scores_mean_df = pd.DataFrame({col: [means[col]] for col in list(single_df.columns)[:-1]})
#         scores_mean_df['disease'] = [disease]
#         # add the means to the results dataframe
#         if results_df is None:
#             results_df = scores_mean_df
#         else:
#             results_df = pd.concat([results_df, scores_mean_df])
#     results_df.to_csv(results_name)
#     # pd.DataFrame(ranked_gene_names).to_csv(results_name + '_ranked_res.csv')

def run_mcl(edgelist, gene_sets, results_name):
    results_df = None
    scores = {'500 count': [], '500 %': [], '100 count': [], '100 %': [], '50 count': [], '50 %': [], '25 count': [],
              '25 %': [], '10 count': [], '10 %': [], 'disease': []}

    disease_gene_sets = {}
    for line in open(gene_sets, 'r'):
        row = line.strip().split('\t')
        # the first item is the disease common name, everything else is genes
        disease_gene_sets[row[0]] = row[1:]
    for disease in disease_gene_sets.keys():
        disease_genes = disease_gene_sets[disease]

        num_edgelist_name = edgelist.split('.')[0] + '_numbered.' + edgelist.split('.')[1]
        mapping = create_numbered_edgelist(edgelist, num_edgelist_name)

        inverse_mapping = {mapping[key]: key for key in mapping.keys()}

        G = nx.read_edgelist(num_edgelist_name)

        matrix = nx.to_scipy_sparse_matrix(G)

        result = mc.run_mcl(matrix)  # run MCL with default parameters
        clusters = mc.get_clusters(result)  # get clusters
        for i in range(50):
            # split the disease set randomly in 2
            nodes = partition(disease_genes, 2)
            seeds, targets = nodes[0], nodes[1]

            # rank genes for based on number of targets in it's cluster
            nodes = []
            target_count = []

            # convert seeds to their numerical equivalent or -1 is they are not in the network
            numbered_seeds = [mapping[x] if x in mapping else -1 for x in seeds if x in mapping]
            numbered_targets = [mapping[x] if x in mapping else -1 for x in targets if x in mapping]

            for c in clusters:
                num_seeds_in_cluster = sum(t in c for t in numbered_seeds)
                for n in c:
                    # make sure that we only include nodes that are genes
                    if 'http' not in inverse_mapping[n]:
                        nodes.append(n)
                        target_count.append(num_seeds_in_cluster)

            df = pd.DataFrame({'nodes': nodes, 'target_count': target_count})
            df = df.sort_values('target_count', ascending=False)

            num_targets_in_network = sum(str(mapping[t]) in list(G.nodes) for t in targets if t in mapping)
            top_xs = [500, 100, 50, 25, 10, 0]
            for i in range(len(top_xs[:-1])):
                print(set(df.iloc[top_xs[i + 1]:top_xs[i], 0]))
                top_overlap = len(list(set(df.iloc[top_xs[i + 1]:top_xs[i], 0]) & set(numbered_targets)))
                scores[str(top_xs[i]) + ' count'].append((top_overlap))
                scores[str(top_xs[i]) + ' %'].append(top_overlap / num_targets_in_network)
            scores['disease'].append(disease)
            single_df = pd.DataFrame(scores)
            # get the mean of each of the iterations for this disease
            means = single_df.mean(axis=0)
            print(means)
            scores_mean_df = pd.DataFrame({col: [means[col]] for col in list(single_df.columns)[:-1]})
            scores_mean_df['disease'] = [disease]
            # add the means to the results dataframe
            if results_df is None:
                results_df = scores_mean_df
            else:
                results_df = pd.concat([results_df, scores_mean_df])
        results_df.to_csv(results_name)


run_mcl(sys.argv[1], sys.argv[2], sys.argv[3])
