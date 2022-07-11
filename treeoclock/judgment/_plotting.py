from treeoclock.clustree.spectral_clustree import _spectral_clustree
import numpy as np
import os
from treeoclock.visualize.plot_coords import plot_coords
from treeoclock.trees.time_trees import TimeTreeSet



# todo include the plot all chains tsne funciton here at some point!


def all_chains_spectral_clustree(cmchain, beta=1, n_clus=2, rf: bool = False, dim: int = 2):
    clustering = cmchain.clustree_all(n_clus=n_clus, beta=beta, rf=rf)
    coords = cmchain.tsne_all(rf=rf, dim=dim)
    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all"
                                                    f"{'_rf' if rf else ''}.png")


def all_chains_added_summaries(cmchain, rf: bool = False, dim: int = 2):
    
    coords = cmchain.tsne_all(rf=rf, dim=dim)

    clustering = np.ones(coords.shape[0])
    # Todo hard coded for one center per chain i.e. half the chains in the given cmchain are centers!
    for i in range(1, int(cmchain.m_MChains / 2)+1):
        clustering[-i] = i+1

    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"summaries{'_rf' if rf else ''}.png", centers=int(cmchain.m_MChains / 2))
