import pandas as pd
from scipy.spatial import distance
import numpy as np

df = pd.read_csv(
    'Embeddings/genes_only_embedding.emb',
    skiprows=1, sep=' ')

# M = np.zeros((df.shape[0], df.shape[0]))
#
# # compute distance over just have the matrix because it is symetric
# for i in range(df.shape[0]):
#     if i % 1000 == 0:
#         print(i)
#     for j in range(df.shape[0] - i):
#         if i == j:
#             continue
#         M[i, j] = distance.euclidean(list(df.iloc[i, 1:]), list(df.iloc[j, 1:]))
#         # make sure the other half of the matrix gets the value too
#         M[j, i] = M[i, j]


def single_row(x):
    i, j = x[0], x[1]
    global df
    a = list(df.iloc[i, 1:])
    b = list(df.iloc[j, 1:])
    return [i,j,distance.euclidean(a,b)]


import multiprocessing as mp

args = []
for i in range(df.shape[0]):
    # A = list(df.iloc[i, 1:])
    for j in range(df.shape[0]-i):
        # B = list(df.iloc[j, 1:])
        args.append([i,j])


pool = mp.Pool(12)
print('pooling')
res = pool.map(single_row,args)
#
# outfile = 'test.tsv'
# start = 0
# with open(outfile, 'w') as out:
#     # compute distance over just have the matrix because it is symetric
#     for i in range(start, start + 1000):
#         print(i)
#         out.write(str(i))
#         A = list(df.iloc[i, 1:])
#         temp = []
#         for j in range(df.shape[0]):
#             if j % 500 == 0:
#                 print('\t' + str(j))
#             temp.append(str(distance.euclidean(A, list(df.iloc[j, 1:]))))
#         out.write("\t".join(temp) + '\n')
