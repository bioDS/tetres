__author__ = "Lars Berling"

from treeoclock.clustree.spectral_clustree import spectral_clustree_dm
from treeoclock.summary.centroid import Centroid

import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import random


def mean_normalized_distance(treeset, summaries, clustering, local=-1):
    """
    Computes the mean normalized distances given a set of trees, set of summary trees and a clustering

    :param trees:
    :type trees: list of list of sets
    :param summaries:
    :type summaries: list of list of sets
    :param clustering: The clustering vecor
    :type clustering: ndarray of int
    :param local: Whether to use the max distance between two trees (local) or the diameter for normalization
    :type local: int
    :return: The mean normalized distance for each cluster
    :rtype: dict[int]->float
    """

    m = len(treeset)  # number of trees
    n = len(treeset[0])  # number of taxa

    # todo optimize this with a given local parameter
    if local != -1:
        # This will normalize with the maximum distance two trees have in the given set
        divisor = local
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


def BIC(treeset, matrix, k, local_normalization, working_folder, add_random=False, chain_id='', _overwrite=False):

    # todo add other clusterings in the future

    # check if working dir exists
    if not os.path.exists(working_folder):
        raise ValueError('Given working dir does not exist')
    # trying to create the cluster folder if not exists
    cluster_folder = os.path.join(working_folder, 'clustering')
    try:
        os.mkdir(cluster_folder)
    except FileExistsError:
        pass

    m = matrix.shape[0]

    summaries = []
    if k is 1:
        clustering = np.zeros(m, dtype=int)

        if _overwrite:
            try:
                os.remove(f'{cluster_folder}/sc-{k}_0.tree')
            except FileNotFoundError:
                pass

        if os.path.exists(f'{cluster_folder}/sc-{k}_0.tree'):
            raise ValueError('NEED to implement read centroid from file, as timetree')
            # todo do i need to check the map for this? probably just assume that it is correct for now
        else:
            # todo call the centroid 5 times and take the best tree out of all of them  -- should be a parameter given
            centroid = Centroid()
            cen, sos = centroid.compute_centroid(treeset)
            with open(f'{cluster_folder}/sc-{k}_0.tree', 'w+') as summary_out:
                # todo write time tree to newick
                # summary_out.write(f'{str(cen.get_newick(f=5))}\n')
                raise ValueError('NEED to implement function, write tree to newick')
        summaries.append(cen)
    else:
        if _overwrite:
            try:
                os.remove(f"{cluster_folder}/sc-{k}-{chain_id}.csv")
            except FileNotFoundError:
                pass

        if os.path.exists(f"{cluster_folder}/sc-{k}-{chain_id}.csv"):
            clustering = np.genfromtxt(f"{cluster_folder}/sc-{k}-{chain_id}.csv", dtype=int)
        else:
            clustering = spectral_clustree_dm(matrix, n_clus=k, beta=1)
            np.save(
                file=f"{cluster_folder}/sc-{k}-{chain_id}.csv",
                # todo maybe add rf option in the future, or other distances
                arr=clustering)

        if add_random:
            random.shuffle(clustering)

        for k_cluster in range(k):
            # todo test this, maybe i need to do something else than just list comprehension here
            cur_treeset = [treeset[index] for index, x in enumerate(clustering) if x == k_cluster]
            # cur_trees = [trees[index] for index, x in enumerate(clustering) if x == cluster]
            if len(cur_treeset) is not 0:  # todo empty clusters????
                if _overwrite:
                    try:
                        os.remove(f'{cluster_folder}/sc-{k}_{k_cluster}.tree')
                    except FileNotFoundError:
                        pass

                if os.path.exists(f'{cluster_folder}/sc-{k}_{k_cluster}.tree'):
                    raise ValueError('Need to read tree at this point!')
                else:
                    # todo change to do best of 5 centroids -- should be a parameter given
                    centroid = Centroid()
                    cen, sos = centroid.compute_centroid(cur_treeset)
                    with open(f'{cluster_folder}/sc-{k}_0.tree', 'w+') as summary_out:
                        # todo write time tree to newick
                        # summary_out.write(f'{str(cen.get_newick(f=5))}\n')
                        raise ValueError('NEED to implement function, write tree to newick ...')
                summaries.append(cen)
            else:
                summaries.append([None])
    local_norm = -1
    if local_normalization:
        local_norm = np.max(matrix)
    mnd_cluster = mean_normalized_distance(treeset, summaries, clustering, local=local_norm)

    # TODO division by 0 is sometimes possible for bic computation,
    #  maybe only if empty cluster exists, therefore delete that bit
    bic = k / 2 * np.log(m) - 2 * m * np.log(1 - np.sum(list(mnd_cluster.values())))
    return bic, mnd_cluster


def plot_bic(treeset, matrix, max_cluster=5, local=False, add_random=False, working_folder="", name=""):
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
    :param local: Whether to use local scale or diameter for normalization of the distance, defaults to False
    :type local: bool
    :param add_random: Add BIC of a randomly shuffled clustering, defaults to False
    :type add_random: bool
    :param ward: Whether to use Ward or spectral clustering, defaults to False
    :type ward: bool
    :return: Saves the BIC cluster evaluation plot
    """

    # mkdir_p(f'MDS_Plots/{tree_name}/Cluster_Summary')
    global max_dist
    bic_values = []
    mnd_list = []
    random_bic = []
    rand_mnd_list = []
    for i in range(max_cluster):
        bic, mnd_cluster = BIC(treeset=treeset,
                           matrix=matrix,
                           k=i+1, local=False,
                           working_folder=working_folder,
                           add_random=False,
                           matrix_is_similarity=False)
        bic_values.append(bic)
        mnd_list.append(mnd_cluster)
        # if add_random:
        #     rand_bic, rand_mnd = BIC(tree_name, i + 1, local, add_random)
        #     random_bic.append(rand_bic)
        #     rand_mnd_list.append(rand_mnd)

    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    df = pd.DataFrame(mnd_list)

    df.plot.bar(stacked=True, ax=ax)

    ax.set_ylabel('Mean normalised distance')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))

    ax2.plot(range(len(bic_values)), bic_values, color='black', label='BIC')
    # if add_random:
    #     ax2.plot(range(len(random_bic)), random_bic, color='blue', label='random BIC')
    
    ax2.set_ylabel('BIC')
    ax2.legend()
    ax2.set_xlabel('Number of clusters k')

    sug_k = bic_values.index(min(bic_values))+1

    plt.suptitle(f"BIC - k = {sug_k} clusters suggested")

    plt.xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')

    # mkdir.mkdir_p("MDS_Plots/BIC_Plots")

    plt.savefig(f"{working_folder}/plots/BIC_{'random_' if add_random else ''}{'local_' if local else ''}{name}.png",
                bbox_inches='tight', dpi=200)


    max_dist = None
    plt.close()

    # Not needed right now
    # todo if add_random is true we should make it a subplot left normal and right the random clustering
    # if add_random:
    #     fig, ax = plt.subplots()
    #     df = pd.DataFrame(rand_mnd_list)
    #     df.plot.bar(stacked=True, ax=ax)
    #     box = ax.get_position()
    #     ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #     ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
    #     plt.xlabel('Number of cluster k')
    #     plt.ylabel('Mean normalised distance')
    #     plt.suptitle('Mean Normalized distance for random clustering')
    #     plt.xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')
    #     plt.savefig(
    #         f"MDS_Plots/BIC_Plots/MND_random_{'local_' if local else ''}{tree_name}.png",
    #         bbox_inches='tight')
    #     plt.clf()
    plt.clf()


def plot_cluster_sizes(tree_name, max_cluster=10):

    """
    Reads all the clusterings from 2 to max_cluster and plots their sizes as a stacked bar plot.

    :param tree_name: Specifies the folder in MDS_Plots/{tree_name}
    :type tree_name: string
    :param max_cluster: The max number for which to plot the sizes, defaults to 10
    :type max_cluster: int
    :return: Saves a cluster sizes plot in MDS_Plots/{tree_name}
    """

    raise ValueError('Not adapted to this package yet ...')

    # cluster_sizes = []
    # for i in range(max_cluster):
    #     if i > 0:
    #         clustering = np.genfromtxt(f'MDS_Plots/{tree_name}/{i + 1}clustering_{tree_name}.csv', dtype=int)
    #     else:
    #         clustering = np.zeros(
    #             len(np.genfromtxt(f'MDS_Plots/{tree_name}/{i + 2}clustering_{tree_name}.csv', dtype=int)),
    #             dtype=int)
    #     cluster_sizes.append(collections.Counter(clustering))
    #
    # df = pd.DataFrame(cluster_sizes)
    #
    # ax = df.plot.bar(stacked=True)
    #
    # plt.xticks(range(max_cluster), labels=[i + 1 for i in range(max_cluster)], rotation='horizontal')
    # plt.ylabel('Number of trees')
    # plt.xlabel('Number of clusters k')
    #
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # # plt.show()
    # plt.savefig(f"MDS_Plots/{tree_name}/Cluster_sizes_{tree_name}.png", bbox_inches='tight')

