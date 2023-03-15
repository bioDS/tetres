from treeoclock.judgment.distance_vs_logparameters import *
from treeoclock.judgment.distance_vs_logparameters import _kc


# todo this testing is nothing that is useful


def test_plot_distance_posterior_comparison(thirty_taxa_MChain):
    plot_distance_posteroir_comparison(thirty_taxa_MChain)
    assert True


def test_plot_sos_log_comparison(thirty_taxa_MChain):
    plot_sos_log_comparison(thirty_taxa_MChain)
    assert True


def test_loglike_along_path(thirty_taxa_MChain):
    plot_loglik_along_path(thirty_taxa_MChain)
    assert True


def test_plot_log_neighbours(thirty_taxa_MChain):
    for dist_type in ["kc"]:  # "rnni", "rf"]:
        plot_log_neighbours(thirty_taxa_MChain, dist_type=dist_type)
    assert True
