import pandas as pd
from scipy.spatial.distance import pdist, squareform
import numpy as np

def hamming(a, b):
    n = a - b
    s = np.count_nonzero(n)
    return s

def L1(a, b):
    n = a - b
    s = sum(np.absolute(n))
    return s
def parse(path):
    cnp = pd.read_csv(path + '/gt.cnp', sep = '\t')
    cnp.drop(cnp.columns[[0,1,2]], axis = 1, inplace = True)
    cnp.columns = cnp.columns.str.replace(' ', '')
    M = cnp.values
    M = M.transpose()
    M[M > 9] = 9
    L = len(M)
    dist = np.zeros((L,L))
    for i in range(0, len(M)):
        for j in range(0, len(M)):
            value = L1(M[i], M[j])
            dist[i][j] = value
            dist[j][i] = value
    #leaf_name = list(cnp.columns.values)
    d = cnp.to_dict('list')
    #d2 = {tuple(v): k for k, v in d.items()}  # exchange keys, values
    #d = {v: list(k) for k, v in d2.items()}   # exchange again
    leaf_name = list(d.keys())
    return path, dist, leaf_name, d
