import numpy as np
from sklearn.cluster import SpectralClustering


def spectral_clustree_dm(matrix, n_clus=2, beta=1):
    similarity = np.exp(-beta * matrix / matrix.std())
    return _spectral_clustree(similarity, n_clus=n_clus)


def _spectral_clustree(matrix, n_clus=2):

    sc = SpectralClustering(n_clusters=n_clus, affinity='precomputed', n_init=100)
    clustering = sc.fit_predict(matrix)

    return np.asarray(clustering, dtype=int)
