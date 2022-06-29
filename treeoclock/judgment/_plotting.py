from treeoclock.clustree.spectral_clustree import _spectral_clustree
import numpy as np
import os
from treeoclock.visualize.plot_coords import plot_coords
from treeoclock.visualize.tsne import _tsne_coords_from_pwd

# todo include the plot all chains tsne funciton here at some point!


def all_chains_spectral_clustree(cmchain, beta=1, n_clus=2):
    # todo only if c is not avaialable, i.e. add saving for clustering
    if not os.path.exists(f"{cmchain.working_dir}/data/{cmchain.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.npy"):
        similarity = cmchain.similarity_matrix_all(beta=beta)
        clustering = _spectral_clustree(similarity, n_clus=n_clus)
        np.save(file=f"{cmchain.working_dir}/data/{cmchain.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.npy",
                arr=clustering)
    else:
        clustering = np.load(f"{cmchain.working_dir}/data/{cmchain.name}_{'' if beta == 1 else f'{beta}_'}clustering_all.npy")

    coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=cmchain.pwd_matrix_all(), dim=2)
    print(kl_divergence)
    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.png")
