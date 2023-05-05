__author__ = "Lars Berling"

import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def mean_normalized_distance(treeset, summaries, clustering, local_norm=False):

    m = len(treeset)  # number of trees
    n = len(treeset[0])  # number of taxa

    # todo optimize this with a given local parameter
    if local_norm:
        # This will normalize with the maximum distance two trees have in the given set
        divisor = local_norm
        # print(f'{divisor} -- max dist\n{((n - 1) * (n - 2) / 2)} -- diameter')
    else:
        # Normalizing with the maximum distance two trees can have, global normalization
        divisor = ((n - 1) * (n - 2) / 2)
    mnd_cluster = {k: 0 for k in range(len(summaries))}
    for i, cluster in enumerate(clustering):
        mnd_cluster[cluster] += treeset[i].fp_distance(summaries[cluster])  # todo this could use multithreading maybe
    for k in mnd_cluster.keys():
        mnd_cluster[k] /= (m * divisor)
    return mnd_cluster


def bic(treeset, clustering, summaries, local_norm):
    # local_norm is either False or the value by which to divide the distances, local = max(distance matrix)

    k = len(summaries)  # number of clusters
    m = len(treeset)  # number of trees in the set

    mnd_cluster = mean_normalized_distance(treeset, summaries, clustering, local_norm=local_norm)
    bic_value = k / 2 * np.log(m) - 2 * m * np.log(1 - np.sum(list(mnd_cluster.values())))

    return bic_value, mnd_cluster


def plot_bic(treeset, clusterings, summaries_dict, max_cluster=5, local_norm=False, add_random=False, file_name = ""):
    # clusterings is a dict k --> k-clustering
    # summaries_dict is a dict k --> list of summaries for k clustering
    # local_norm is either False or the value used to normalize distances

    bic_values = []
    mnd_list = []
    random_bic = []
    rand_mnd_list = []
    for i in range(max_cluster):
        bic_value, mnd_cluster = bic(treeset=treeset,
                                       clustering=clusterings[i+1],
                                       summaries=summaries_dict[i+1],
                                       local_norm=local_norm)

        bic_values.append(bic_value)
        mnd_list.append(mnd_cluster)
        if add_random:
            rshuf_cluster = clusterings[i+1].copy()
            random.shuffle(rshuf_cluster)
            rand_bic, rand_mnd = bic(treeset=treeset,
                                       clustering=rshuf_cluster,
                                       summaries=summaries_dict[i+1],
                                       local_norm=local_norm)
            random_bic.append(rand_bic)
            rand_mnd_list.append(rand_mnd)

    df = pd.DataFrame(mnd_list)  # creating mean normalized distances dataframe

    if not add_random:
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        df.plot.bar(stacked=True, ax=ax)
        ax.set_ylabel('Mean normalised distance')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
        ax2.plot(range(len(bic_values)), bic_values, color='black', label='BIC')
        ax2.set_ylabel('BIC')
        ax2.set_xlabel('Number of clusters k')
        mnd_constraint = np.argmin([sum(c.values()) for c in mnd_list]) + 1
        sug_k = np.argmin(bic_values[0:mnd_constraint])+1
        plt.suptitle(
            f"BIC - k = {sug_k} clusters suggested {'(MND constraint)' if mnd_constraint != max_cluster else ''}")
        plt.xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')

    else:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 7))  # todo figure out a good size
        ax2 = ax[0].twinx()
        df.plot.bar(stacked=True, ax=ax[0])
        ax[0].set_ylabel('Mean normalised distance')
        box = ax[0].get_position()
        ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax[0].legend(loc='center left', bbox_to_anchor=(1.15, 0.5))
        ax2.plot(range(len(bic_values)), bic_values, color='black', label='BIC')
        ax2.set_ylabel('BIC')
        ax2.set_xlabel('Number of clusters k')
        mnd_constraint = np.argmin([sum(c.values()) for c in mnd_list]) + 1
        sug_k = np.argmin(bic_values[0:mnd_constraint]) + 1
        plt.suptitle(
            f"BIC - k = {sug_k} clusters suggested {'(MND constraint)' if mnd_constraint != max_cluster else ''}")
        ax[0].set_xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')

        # adding the random plot to the second subplot
        df_rand = pd.DataFrame(rand_mnd_list)
        ax3 = ax[1].twinx()
        df_rand.plot.bar(stacked=True, ax=ax[1])
        box = ax[1].get_position()
        ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax3.plot(range(len(random_bic)), random_bic, color='blue', label='random BIC')
        ax3.set_ylabel("BIC")
        ax[1].set_ylabel('Mean normalised distance')
        ax[1].legend(loc='center left', bbox_to_anchor=(1.15, 0.5))
        ax3.set_xlabel('Number of cluster k')
        ax[1].title.set_text("Random clustering")
        ax[1].set_xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')

    # File format of the plot is inferred from the name, i.e. .pdf/.eps/.png
    # Todo maybe make the plot parameters an argument with kwargs in the future
    plt.tight_layout()
    if file_name != "":
        plt.savefig(file_name, bbox_inches='tight')
        plt.clf()
        plt.close()
    else:
        plt.show()
