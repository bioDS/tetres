from tetres.clustree.bic import bic
from tetres.trees.time_trees import TimeTreeSet
import numpy as np
import os


def bic_test(multichain, i, j):
    # sorting i and j
    if i > j:
        i, j = j, i
    # Building the combined matrix from the pwd matrix function
    row1 = multichain.pwd_matrix(i)
    row1 = np.concatenate((row1, multichain.pwd_matrix(i, j)), axis=1)
    row2 = np.zeros((len(multichain[i].trees), len(multichain[j].trees)))
    row2 = np.concatenate((row2, multichain.pwd_matrix(j)), axis=1)
    matrix = np.concatenate((row1, row2), axis=0)

    full_trees = TimeTreeSet()
    full_trees.trees = multichain[i].trees.copy()
    full_trees.trees.extend(multichain[j].trees.copy())
    full_trees.map = multichain[i].trees.map

    # calculating the pooled BIC for all trees in one cluster
    bic_pooled = bic(treeset=full_trees,
                     k=1,
                     chain_id=f"{multichain.name}_bic_test({i},{j})",
                     matrix=matrix,
                     working_folder=multichain.working_dir,
                     local_norm=False)
    # temporary saving of a 2 clustering that is appropriate for this test, i.e. each chain is 1 cluster
    cluster_folder = os.path.join(multichain.working_dir, 'clustering')

    clustering = np.concatenate((np.zeros(len(multichain[i].trees), dtype=int),np.ones(len(multichain[j].trees), dtype=int)), axis=0)
    np.save(file=f"{cluster_folder}/sc-2-{multichain.name}_bic_test({i},{j})", arr=clustering)

    bic_clustered = bic(treeset=full_trees,
                        k=2,
                        chain_id=f"{multichain.name}_bic_test({i},{j})",
                        matrix=matrix,
                        working_folder=multichain.working_dir,
                        local_norm=False)
    # todo maybe check something with the MND
    #  Currently only the BICs are returned
    return bic_pooled[0], bic_clustered[0]