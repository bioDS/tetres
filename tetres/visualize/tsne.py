import numpy as np
from sklearn import manifold
import os


# todo this should go into the clustree folder/subpackage

def _tsne_coords_from_pwd(pwd_matrix, dim=2):
    internal_matrix = pwd_matrix + pwd_matrix.transpose()
    mds = manifold.TSNE(n_components=dim,
                        metric="precomputed",
                        learning_rate='auto',
                        init="random",
                        perplexity=10)
    coords = mds.fit_transform(internal_matrix)
    return coords, mds.kl_divergence_


def tsne_coords_from_mchain(mchain, dim=2):
    # computing or loading the tsne coordinates
    # todo maybe change this to do a recomputation if wanted and overwriting of the old if better fit

    # todo this should just be a funciton that computes tsne, all the questions should be part of the mchain function

    if not os.path.exists(f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy"):
        pwd_matrix = mchain.pwd_matrix()
        coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=pwd_matrix, dim=dim)
        # todo do something with the kl_divergence value
        np.save(file=f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy", arr=coords)
        return coords
    else:
        return np.load(f"{mchain.working_dir}/data/{mchain.name}_tsne_{dim}d_coords.npy")
