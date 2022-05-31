from treeoclock.visualize.tsne import tsne_coords_from_mchain
from treeoclock.visualize.plot_coords import plot_coords, plot_density_over_coordinates

import numpy as np


def test_plot_coords_2d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=2)
    plot_coords(coords, colors=np.asarray(ten_taxa_cMChain[1].log_data["likelihood"]))


def test_plot_coords_3d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=3)
    plot_coords(coords, colors=np.asarray(ten_taxa_cMChain[1].log_data["posterior"]))


def test_plot_density(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[1], dim=2)
    plot_density_over_coordinates(coords, np.asarray(ten_taxa_cMChain[1].log_data["posterior"]),
                                  filename=f"{ten_taxa_cMChain.working_dir}/testing/")