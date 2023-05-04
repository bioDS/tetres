__author__ = "Lars Berling"


from tetres.trees.time_trees import TimeTreeSet

import numpy as np
import matplotlib.pyplot as plt
import os.path
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
    mnd_cluster = {k: 0 for k in range(max(clustering) + 1)}
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


def plot_bic(treeset, matrix, max_cluster=5, local_norm=False, add_random=False, working_folder="", name="", _overwrite=False):
    """
    Saves the BIC plot for cluster evaluation

    Computes the Bayesian inference criterion for the clusterings in the rango of 1:max_cluster.
    The parameters specify: Whether the mean normalized distance should use the diameter or the max distance
    of two trees for normalization, whether to add a randomly shuffled clustering
    and which clustering to use (Ward or spectral clustering).

    :param tree_name: Specifies the folder MDS_Plots/{tree_name}
    :type tree_name: string
    :param max_cluster: The max number for which to plot/compute the BIC, defaults to 10
    :type max_cluster: int
    :param local_norm: Whether to use local scale or diameter for normalization of the distance, defaults to False
    :type local_norm: bool
    :param add_random: Add BIC of a randomly shuffled clustering, defaults to False
    :type add_random: bool
    :return: Saves the BIC cluster evaluation plot
    """
    if _overwrite:
        try:
            os.remove(f"{working_folder}/plots/BIC_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png")
        except FileNotFoundError:
            pass
    if os.path.exists(f"{working_folder}/plots/BIC_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png"):
        raise FileExistsError("Plot already exists -- pass _overwrite=True if you wish to recompute.")

    bic_values = []
    mnd_list = []
    random_bic = []
    rand_mnd_list = []
    for i in range(max_cluster):
        bic_value, mnd_cluster = bic(treeset=treeset,
                               matrix=matrix,
                               k=i+1, local_norm=local_norm,
                               working_folder=working_folder,
                               _overwrite=_overwrite,
                               chain_id=name
                               )
        bic_values.append(bic_value)
        mnd_list.append(mnd_cluster)
        if add_random:
            rand_bic, rand_mnd = bic(treeset=treeset,
                                     matrix=matrix,
                                     k=i+1, local_norm=local_norm,
                                     working_folder=working_folder,
                                     _overwrite=_overwrite,
                                     chain_id=name,
                                     random_shuffle=True)
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
        sug_k = bic_values.index(min(bic_values))+1
        plt.suptitle(f"BIC - k = {sug_k} clusters suggested")
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
        sug_k = bic_values.index(min(bic_values)) + 1
        plt.suptitle(f"BIC - k = {sug_k} clusters suggested")
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

    # todo change to eps, or make the option available at some point
    plt.tight_layout()
    plt.savefig(f"{working_folder}/plots/BIC_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png",
                bbox_inches='tight', dpi=200)
    plt.clf()
    plt.close()
