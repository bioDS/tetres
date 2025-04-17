from tetres.visualize.plot_config import PlotOptions
from tetres.visualize.plot_coords import plot_coords

import numpy as np


# todo colors as keyword for log_data and generally the function should be called from the class
#  instead of directly calling the function
def test_plot_coords_2d(ten_taxa_multichain):
    coords = ten_taxa_multichain.get_mds_coords(target=1, _overwrite=True)
    plot_options = PlotOptions(
        colors=np.asarray(ten_taxa_multichain[1].log_data["likelihood"])
    )
    plot_coords(coords, dim=2, options=plot_options)


def test_plot_coords_2d_all(ten_taxa_multichain):
    coords = ten_taxa_multichain.get_mds_coords(target="all", _overwrite=True)
    plot_options = PlotOptions(
        colors=np.concatenate(
            [[i] * len(chain) for i, chain in enumerate(ten_taxa_multichain.MChain_list)]),
        label=np.concatenate(
            [[ten_taxa_multichain[i].name] * len(chain)
             for i, chain in enumerate(ten_taxa_multichain.MChain_list)]),
        title="TSNE Plot"
    )
    plot_coords(coords, dim=2, options=plot_options)

