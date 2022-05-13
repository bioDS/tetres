# Analysis of tree sets that are interpreted as a true posteior or an estimate of that

from treeoclock.judgment.mchain import MChain

import itertools
import operator
import os
from ctypes import CDLL, POINTER, c_long
from treeoclock.trees._ctrees import TREE_LIST, TREE
import numpy as np


lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


def calc_pw_distances(trees):

    # todo add a n_cores parameter to c function?

    n = len(trees)
    lib.pairwise_distance.argtypes = [POINTER(TREE_LIST), c_long, POINTER((c_long * n) * n)]
    distances = np.zeros((n, n), dtype=c_long)
    ctreelist = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees]))
    lib.pairwise_distance(ctreelist, n, distances.ctypes.data_as(POINTER((c_long * n) * n)))
    return distances


def nbr_unique_trees(mchain: MChain, dm_name=""):

    # todo add the pairwise distance matrix as part of the MChain object, makes a lot of methods faster if recomputed or in general faster i believe!

    if not os.path.exists(f"{mchain.working_dir}/{dm_name}.csv"):
        pd = calc_pw_distances(mchain.trees)
        if dm_name:
            # todo this needs to be more efficient!!!
            np.savetxt(f"{mchain.working_dir}/{dm_name}.csv", pd, delimiter=',', fmt='%i')
    else:
        # todo this should also be more efficient!
        pd = np.genfromtxt(f"{mchain.working_dir}/{dm_name}.csv", delimiter=',', dtype=int)
        # pd.astype(int, inplace=True)
    
    unique = 0
    if not os.path.exists(f"{mchain.working_dir}/{dm_name}_uniques.csv"):
        count = 0
        remove = set()
        for i in range(0, len(mchain.trees) - 1):
            zero = False
            for j in range(i + 1, len(mchain.trees)):
                if pd[i, j] == 0:
                    zero = True
                    remove.add(j)
            if zero:
                count += 1
        unique = np.delete(pd, list(remove), 0)  # todo there is a problem with this deletion part i think
        np.savetxt(f"{mchain.working_dir}/{dm_name}_uniques.csv", unique, delimiter=',', fmt='%i')
    else:
        unique = np.genfromtxt(f"{mchain.working_dir}/{dm_name}_uniques.csv", delimiter=',', dtype=int)
        # unique.astype(int, inplace=True)
    
    print(len(unique))
    # todo check if the set is connected in RNNI based on BFS over the matrix
    #  put node 0 in set and all nodes at distance 1,

    return 0
    
    
