import numpy as np
from sklearn.cluster import SpectralClustering

###
# Simply calculate the spectral clustering from a distance matrix,
# will always calculate the similarity matrix
###


def spectral_clustree_dm(matrix, n_clus=2, beta=1):
    if np.allclose(matrix, np.triu(matrix)):
        # matrix is upper triangular, has to be changed to a full matrix for this implementation
        matrix += np.transpose(matrix)
    similarity = np.exp(-beta * matrix / matrix.std())
    return _spectral_clustree(similarity, n_clus=n_clus)


def _spectral_clustree(matrix, n_clus=2):
    if np.allclose(matrix, np.triu(matrix)):
        # matrix is upper triangular, has to be changed to a full matrix for this implementation
        matrix += np.transpose(matrix)
    sc = SpectralClustering(n_clusters=n_clus, affinity='precomputed', n_init=100)
    clustering = sc.fit_predict(matrix)
    return np.asarray(clustering, dtype=int)
