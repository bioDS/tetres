from tetres.visualize.tsne import tsne_coords_from_mchain
from tetres.visualize.plot_coords import plot_coords, plot_density_over_coordinates, plot_all_chains_tsne

import numpy as np


# todo colors as keyword for log_data and generally the function should be called from the class
#  instead of directly calling the function
def test_plot_coords_2d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=2)
    plot_coords(coords, colors=np.asarray(ten_taxa_cMChain[1].log_data["likelihood"]))


def test_plot_coords_3d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=3)
    plot_coords(coords, colors=np.asarray(ten_taxa_cMChain[1].log_data["posterior"]))


def test_plot_coords_cluster_2d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=2)
    plot_coords(coords, colors=ten_taxa_cMChain[1].spectral_clustree(), colorbar=True, centers=["C1", "C2"], scale_centers=3)


def test_plot_coords_cluster_3d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=3)
    plot_coords(coords, colors=ten_taxa_cMChain[1].spectral_clustree(), colorbar=True, centers=["C1", "C3"], scale_centers=3)


def test_plot_density(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=2)
    plot_density_over_coordinates(coords, np.asarray(ten_taxa_cMChain[1].log_data["posterior"]),
                                  filename=f"{ten_taxa_cMChain.working_dir}/plots")


def test_plot_all_chains_tsne(ten_taxa_cMChain):
    plot_all_chains_tsne(ten_taxa_cMChain)
    assert True
