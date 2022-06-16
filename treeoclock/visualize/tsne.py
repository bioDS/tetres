import numpy as np
from sklearn import manifold
import os


def _tsne_coords_from_pwd(pwd_matrix, dim=2):
    internal_matrix = pwd_matrix + pwd_matrix.transpose()
    mds = manifold.TSNE(n_components=dim, metric="precomputed", learning_rate='auto')
    coords = mds.fit_transform(internal_matrix)
    return coords


def tsne_coords_from_mchain(mchain, dim=2):
    # computing or loading the tsne coordinates
    if not os.path.exists(f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy"):
        pwd_matrix = mchain.pwd_matrix()
        coords = _tsne_coords_from_pwd(pwd_matrix=pwd_matrix, dim=dim)
        np.save(file=f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy", arr=coords)
        return coords
    else:
        return np.load(f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy")
