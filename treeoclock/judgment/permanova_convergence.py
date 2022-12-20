import gc
import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

global _gr_boundary
_gr_boundary = 0.02

from skbio.stats.distance import permanova
from skbio import DistanceMatrix


def permanova_cut(cmchain, i, j, smoothing=0.9):

    global _gr_boundary

    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    for cur_sample in range(10,dm_i.shape[0], 50):
        slide_start = int(np.rint(cur_sample * (1 - smoothing)))  # TODO smoothing is set to 0.9 for now, only burnin detection for this case!

        cur_block = np.block([[dm_i[slide_start:cur_sample+1, slide_start:cur_sample+1],
                               dm_ij[slide_start:cur_sample+1, slide_start:cur_sample+1]],
                              [dm_ji[slide_start:cur_sample+1, slide_start:cur_sample+1],
                               dm_j[slide_start:cur_sample+1, slide_start:cur_sample+1]]])
        cur_groups = [0 for _ in range(slide_start, cur_sample+1)] + [1 for _ in range(slide_start, cur_sample+1)]
        cur_block = DistanceMatrix(cur_block + cur_block.transpose())
        pm = permanova(cur_block, cur_groups, permutations=1000)

        if pm["p-value"] < 0.05:
            return slide_start, cur_sample

    # No cutoff found
    return -1, -1


def permanova_pvalue_plot(cmchain, i, j, smoothing_l=[0.9]):

    global _gr_boundary

    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    data = []

    for cur_sample in range(10,dm_i.shape[0], 50):
        for smoothing in smoothing_l:
            slide_start = int(np.rint(cur_sample * (1 - smoothing)))  # TODO smoothing is set to 0.9 for now, only burnin detection for this case!

            cur_block = np.block([[dm_i[slide_start:cur_sample+1, slide_start:cur_sample+1],
                                   dm_ij[slide_start:cur_sample+1, slide_start:cur_sample+1]],
                                  [dm_ji[slide_start:cur_sample+1, slide_start:cur_sample+1],
                                   dm_j[slide_start:cur_sample+1, slide_start:cur_sample+1]]])
            cur_groups = [0 for _ in range(slide_start, cur_sample+1)] + [1 for _ in range(slide_start, cur_sample+1)]
            cur_block = DistanceMatrix(cur_block + cur_block.transpose())
            pm = permanova(cur_block, cur_groups, permutations=1000)

            data.append([pm["p-value"], smoothing, cur_sample])
    data = pd.DataFrame(data, columns=["p-value", "smoothing", "sample"])
    sns.swarmplot(data=data, x="sample", y="p-value", hue="smoothing")
    plt.xticks(rotation=270)
    plt.show()
    return 0
