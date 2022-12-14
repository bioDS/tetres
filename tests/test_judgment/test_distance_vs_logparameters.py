from treeoclock.judgment.distance_vs_logparameters import *


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
    plot_log_neighbours(thirty_taxa_MChain, rf=True)
    plot_log_neighbours(thirty_taxa_MChain, rf=False)
    assert True
