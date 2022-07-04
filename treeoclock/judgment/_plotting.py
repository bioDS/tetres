from treeoclock.clustree.spectral_clustree import _spectral_clustree
import numpy as np
import os
from treeoclock.visualize.plot_coords import plot_coords


# todo include the plot all chains tsne funciton here at some point!


def all_chains_spectral_clustree(cmchain, beta=1, n_clus=2, rf: bool = False, dim: int = 2):
    clustering = cmchain.clustree_all(n_clus=n_clus, beta=beta, rf=rf)
    coords = cmchain.tsne_all(rf=rf, dim=dim)
    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all"
                                                    f"{'_rf' if rf else ''}.png")
