# Analysis of tree sets that are interpreted as a true posteior or an estimate of that

from treeoclock.judgment.mchain import MChain
from treeoclock.visualize.graphs import visualize_matrix

import os

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.sparse.csgraph import connected_components
from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances

# todo testing required!
def nbr_unique_trees(mchain: MChain, dm_name=""):
    # todo add the pairwise distance matrix as part of the MChain object, makes a lot of methods faster if recomputed or in general faster i believe!

    if not os.path.exists(f"{mchain.working_dir}/{dm_name}_uniques.csv.gz"):
        if not os.path.exists(f"{mchain.working_dir}/{dm_name}.csv.gz"):
            pd = calc_pw_distances(mchain.trees)
            if dm_name:
                # todo this needs to be more efficient!!!
                np.savetxt(fname=f"{mchain.working_dir}/{dm_name}.csv.gz", X=pd, delimiter=',', fmt='%i')
        else:
            # todo this should also be more efficient!
            pd = np.genfromtxt(fname=f"{mchain.working_dir}/{dm_name}.csv.gz", delimiter=',', dtype=int)
            # pd.astype(int, inplace=True)

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
        unique = np.delete(pd, list(remove), 0)
        unique = np.delete(unique, list(remove), 1)

        np.savetxt(fname=f"{mchain.working_dir}/{dm_name}_uniques.csv.gz", X=unique, delimiter=',', fmt='%i')
    else:
        unique = np.genfromtxt(fname=f"{mchain.working_dir}/{dm_name}_uniques.csv.gz", delimiter=',', dtype=int)
        # unique.astype(int, inplace=True)

    # todo just double check at some point that unique contains what i want it to be!
    n_taxa = len(mchain.trees[0])
    diam = (n_taxa * (n_taxa - 1)) /2

    distance_distribution_matrix(unique, diam)
    indexes = extract_largest_cc_indices(unique.copy())
    distance_distribution_matrix(unique, diam, indexes)

    # visualize_matrix(unique, f"{mchain.working_dir}/largest_cc.dot", f"{mchain.working_dir}/largest_cc.png")

    # bfs = bfs_distance_matrix(d_matrix=unique)

    return 0

from collections import Counter

# todo testing required!
def extract_largest_cc_indices(d_matrix):
    d_matrix[d_matrix > 1] = 0
    ccs = connected_components(d_matrix)
    # print(len(ccs[1]))
    # print(len(set(ccs[1])))
    ccs_count = Counter(ccs[1])
    # print(ccs_count)
    val = max(ccs_count, key=lambda key: ccs_count[key])
    # print(val)
    ind = np.where(ccs[1] == val)[0]
    # print(ind)
    # print(len(ind))
    return ind

# todo testing required!
def distance_distribution_matrix(d_matrix, diam, indexes=None):
    if indexes is not None:
        d_matrix = d_matrix[np.ix_(indexes, indexes)]
    values = d_matrix[d_matrix > 0]
    sns.histplot(values, stat="density", discrete=True)
    plt.axvline(x=diam, color="red", linestyle="--")
    plt.axvline(x=np.mean(values), color="purple", linestyle="-.")
    plt.show()
    plt.clf()


