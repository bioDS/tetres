__author__ = "Lars Berling"

from tetres.clustree.spectral_clustree import spectral_clustree_dm
import numpy as np
import matplotlib.pyplot as plt
import random
import os
from sklearn.metrics import silhouette_score as sksilhouette


def silhouette_score(matrix, k, working_folder, local_norm=False, random_shuffle=False, chain_id='myChain', _overwrite=False):

    if not os.path.exists(working_folder):
        raise ValueError('Given working dir does not exist')
    # trying to create the cluster folder if not exists
    cluster_folder = os.path.join(working_folder, 'clustering')
    try:
        os.mkdir(cluster_folder)
    except FileExistsError:
        pass

    if k == 1:
        raise ValueError("Silhouette for a single cluster not implemented!")

    # todo the clustering should not be calculated here, but rather a parameter, therefore this will be removed soon
    if _overwrite:
        try:
            os.remove(f"{cluster_folder}/sc-{k}-{chain_id}.npy")
        except FileNotFoundError:
            pass
    if os.path.exists(f"{cluster_folder}/sc-{k}-{chain_id}.npy"):
        clustering = np.load(f"{cluster_folder}/sc-{k}-{chain_id}.npy")
    else:
        clustering = spectral_clustree_dm(matrix, n_clus=k, beta=1)
        np.save(file=f"{cluster_folder}/sc-{k}-{chain_id}", arr=clustering)
    if random_shuffle:
        random.shuffle(clustering)
    if np.allclose(matrix, np.triu(matrix)):
        # matrix is upper triangular, has to be changed to a full matrix for this implementation
        matrix += np.transpose(matrix)
    if not local_norm:
        return sksilhouette(metric="precomputed", X=np.asarray(matrix, dtype=float), labels=np.asarray(clustering, dtype=float))

    # Implementation to divide with a local scale instead of the max(a, b) at every steps
    m = matrix.shape[0]  # number of trees
    cluster_indices = []
    for i in range(k):
        # creating list of list of indeces for each cluster
        cluster_indices.append(np.where(clustering == i)[0])

    # initializing everything with 0
    a = np.zeros(m, float)
    b = np.zeros(m, float)
    s = np.zeros(m, float)
    for i in range(m):
        i_clus = clustering[i]

        # Computation of a starts
        intra_clust_dist = 0
        for j in cluster_indices[i_clus]:
            # get matrix entry
            if not i == j:
                intra_clust_dist += matrix[i][j]
        if len(cluster_indices[i_clus]) > 1:
            a[i] = (1 / (len(cluster_indices[i_clus]) - 1)) * intra_clust_dist
        else:
            a[i] = 0
        # Computation of b starts
        cur_b = -1
        for other_cluster in range(k):
            if not other_cluster == i_clus and len(cluster_indices[other_cluster]) > 0:
                # Loop over all other clusters that the one i is in
                inter_clus_dist = float(0)
                for j in cluster_indices[other_cluster]:
                    inter_clus_dist += matrix[i][j]
                inter_clus_dist = (1 / len(cluster_indices[other_cluster])) * inter_clus_dist
                if inter_clus_dist < cur_b or cur_b == -1:
                    cur_b = inter_clus_dist
        if cur_b == -1:
            # Only one cluster is not empty therefore b is 0
            cur_b = 0
        b[i] = cur_b

        if len(cluster_indices[i_clus]) == 1:
            s[i] = 0
        else:
            # s[i] = (b[i] - a[i]) / np.maximum(a[i], b[i])  # The default silhouette score, identical to sklearn implementation
            s[i] = (b[i] - a[i])

    divisor = 0
    for i in range(m):
        diff = np.abs(a[i] - b[i])
        if diff > divisor:
            divisor = diff
    s /= divisor
    return np.mean(s)


def plot_silhouette(matrix, max_cluster=5, local_norm=False, add_random=False, working_folder="", _overwrite=False, name="MyChain"):
    if os.path.exists(f"{working_folder}/plots/Silhouette_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png"):
        if _overwrite:
            os.remove(f"{working_folder}/plots/Silhouette_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png")
        else:
            raise FileExistsError("Use _overwrite=True to overwrite the already existing plot.")

    sil_scores = []
    random_sil_scores = []
    for i in range(2, max_cluster + 1):
        s = silhouette_score(matrix=matrix,
                           k=i, local_norm=True,
                           working_folder=working_folder,
                           random_shuffle=False)
        sil_scores.append(s)
        if add_random:
            rs = silhouette_score(matrix=matrix,
                           k=i, local_norm=True,
                           working_folder=working_folder,
                           random_shuffle=True)
            random_sil_scores.append(np.mean(rs))

    plt.plot(range(2, max_cluster + 1), sil_scores, color='black', label='Sil. Coeff.')
    if add_random:
        plt.plot(range(2, max_cluster + 1), random_sil_scores, color='blue', label='random Sil. Coeff.')

    sug_k = sil_scores.index(max(sil_scores)) + 2

    plt.xlabel('Number of Cluster')
    plt.ylabel('Silhouette Coefficient')
    plt.suptitle(f'Silhouette Coefficient, k = {sug_k} clusters suggested')

    plt.xticks(range(2, max_cluster + 1), labels=[i for i in range(2, max_cluster + 1)], rotation='horizontal')

    plt.legend(loc="upper center", fancybox=True, bbox_to_anchor=(0.5, 1.15), shadow=True, ncol=2)
    plt.tight_layout()

    plt.savefig(f"{working_folder}/plots/Silhouette_{'random_' if add_random else ''}{'local_' if local_norm else ''}{name}.png", bbox_inches='tight', dpi=200)
    plt.clf()
    plt.close()