from treeoclock.clustree.spectral_clustree import _spectral_clustree
import numpy as np
import os
from treeoclock.visualize.plot_coords import plot_coords
from treeoclock.visualize.tsne import _tsne_coords_from_pwd

# todo include the plot all chains tsne funciton here at some point!


def all_chains_spectral_clustree(cmchain, beta=1, n_clus=2, rf: bool = False):
    clustering = cmchain.clustree_all(n_clus=n_clus, beta=beta, rf=rf)

    # todo the coords should be saved, together with the kl_divergence
    coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=cmchain.pwd_matrix_all(rf=rf), dim=2)

    print(kl_divergence)
    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all"
                                                    f"{'_rf' if rf else ''}.png")
